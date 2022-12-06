"""
Salt column test
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest

from credo.jobrunner import SimpleJobRunner
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from credo.systest import FieldWithinTolTC, OneDSolutionWithinTolTC

from mulgrids import mulgrid
from numpy import polyval

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'

import numpy as np
from docutils.core import publish_file

parser = argparse.ArgumentParser()
parser.add_argument("-np", type = int, default = 1, help = "number of processes")
parser.add_argument("-d", "--docker", action = "store_true",
                    help = "run via Docker (waiwera-dkr)")
args = parser.parse_args()
mpi = args.np > 1 and not args.docker
simulator = 'waiwera-dkr -np %d' % args.np if args.docker else 'waiwera'

model_name = 'salt_column'

def AUTOUGH2_solid_saturation(mResult, index):
    """Returns solid saturation from AUTOUGH2 results."""
    sl = mResult.getFieldAtOutputIndex('Liquid satur', index)
    sv = mResult.getFieldAtOutputIndex('Gas saturati', index)
    return 1. - sl - sv

AUTOUGH2_FIELDMAP = {
    'Pressure': 'Pressure',
    'Temperature': 'Temperature',
    'Vapour saturation': 'Gas saturati',
    'Liquid saturation': 'Liquid satur',
    'Liquid salt mass fraction': 'NaCl liquid',
    'Solid saturation': AUTOUGH2_solid_saturation
}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Liquid saturation': 'fluid_liquid_saturation',
    'Liquid salt mass fraction': 'fluid_liquid_salt_mass_fraction',
    'Solid saturation': 'fluid_solid_saturation'
}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
map_out_atm = list(range(geo.num_atmosphere_blocks, geo.num_blocks))

test_fields = ["Pressure", "Temperature", "Liquid saturation", "Vapour saturation",
               "Solid saturation", "Liquid salt mass fraction"]
plot_fields = test_fields
field_scale = {"Pressure": 1.e5}
field_unit = {"Pressure": "bar", "Temperature": "$^{\circ}$C"}

sat_tol = 0.05
field_tols = {"Pressure": 0.01, "Temperature": 0.02, "Liquid saturation": sat_tol,
             "Vapour saturation": sat_tol, "Solid saturation": sat_tol,
             "Liquid salt mass fraction": 0.01}

salt_column_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
salt_column_test.description = """Vertical column salt test. The AUTOUGH2 EWASG EOS uses different brine thermodynamics, so an exact match is not expected."""

run_index = 0
run_name = 'run'
run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                            fieldname_map = WAIWERA_FIELDMAP,
                            simulator = simulator,
                            basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = args.np
salt_column_test.mSuite.addRun(model_run, run_name)

salt_column_test.setupEmptyTestCompsList()

AUTOUGH2_result = {}

run_base_name = model_name
results_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result[run_name] = T2ModelResult("AUTOUGH2", results_filename,
                                          geo_filename = t2geo_filename,
                                          fieldname_map = AUTOUGH2_FIELDMAP,
                                          ordering_map = map_out_atm)
salt_column_test.addTestComp(run_index, "AUTOUGH2",
                  FieldWithinTolTC(fieldsToTest = test_fields,
                                   fieldTols = field_tols,
                                   expected = AUTOUGH2_result[run_name],
                                   testOutputIndex = -1))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = salt_column_test.runTest(jrunner, createReports = True)

z = [lay.centre for lay in geo.layerlist[1:]]
for field_name in plot_fields:
    scale = field_scale[field_name] if field_name in field_scale else 1.
    unit = field_unit[field_name] if field_name in field_scale else ''
    result = salt_column_test.mSuite.resultsList[run_index]
    var = result.getFieldAtOutputIndex(field_name, -1) / scale
    plt.plot(var, z, '-', label = 'Waiwera', zorder = 2)
    var = AUTOUGH2_result[run_name].getFieldAtOutputIndex(field_name, -1) / scale
    plt.plot(var, z, 's', label = 'AUTOUGH2', zorder = 1)
    plt.xlabel(field_name + ' (' + unit + ')')
    plt.ylabel('elevation (m)')
    plt.title(field_name)
    img_filename_base = '_'.join((model_name, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(salt_column_test.mSuite.runs[run_index].basePath,
                                salt_column_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.legend(loc = 'best')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    salt_column_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], salt_column_test.outputPathBase):
    report_filename = os.path.join(salt_column_test.outputPathBase,
                     "%s-report.%s" % (salt_column_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(salt_column_test, mResults, rGen, report_filename)
    html_filename = os.path.join(salt_column_test.outputPathBase,
                     "%s-report.%s" % (salt_column_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
