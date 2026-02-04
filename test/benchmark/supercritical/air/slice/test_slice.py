"""
Supercritical natural-state air/water 2D slice problem
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
from credo.modelresult import DigitisedOneDFieldResult

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

model_name = 'slice'

AUTOUGH2_FIELDMAP = {
    'Pressure': 'Pressure',
    'Temperature': 'Temperature'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
map_out_atm = list(range(geo.num_atmosphere_blocks, geo.num_blocks))

test_fields = ['Pressure', 'Temperature']
field_scale = {'Pressure': 1e6, 'Temperature': 1.}
field_unit = {'Pressure': 'MPa', 'Temperature': '$^{\circ}$C'}
field_tols = {"Pressure": 0.03, "Temperature": 0.01}

slice_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
slice_test.description = """Supercritical natural-state air/water 2D slice problem
"""

run_index = 0
run_name = 'run'
run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                            fieldname_map = WAIWERA_FIELDMAP,
                            simulator = simulator,
                            basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = args.np
slice_test.mSuite.addRun(model_run, run_name)

slice_test.setupEmptyTestCompsList()

AUTOUGH2_result = {}

run_base_name = model_name
results_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result[run_name] = T2ModelResult("AUTOUGH2", results_filename,
                                          geo_filename = t2geo_filename,
                                          fieldname_map = AUTOUGH2_FIELDMAP,
                                          ordering_map = map_out_atm)
slice_test.addTestComp(run_index, "AUTOUGH2",
                       FieldWithinTolTC(fieldsToTest = test_fields,
                                        fieldTols = field_tols,
                                        expected = AUTOUGH2_result[run_name],
                                        testOutputIndex = -1))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = slice_test.runTest(jrunner, createReports = True)

z = [lay.centre for lay in geo.layerlist[1:]]
blknames = [geo.block_name(lay.name, geo.columnlist[0].name) for lay in geo.layerlist[1:]]
blkindices = np.array([geo.block_name_index[blkname] for blkname in blknames],
                      dtype = int) - geo.num_atmosphere_blocks

for field_name in test_fields:
    scale = field_scale[field_name] if field_name in field_scale else 1.
    unit = field_unit[field_name] if field_name in field_scale else ''
    result = slice_test.mSuite.resultsList[run_index]
    var = result.getFieldAtOutputIndex(field_name, -1) / scale
    plt.plot(var[blkindices], z, '-', label = 'Waiwera')
    var = AUTOUGH2_result[run_name].getFieldAtOutputIndex(field_name, -1) / scale
    plt.plot(var[blkindices], z, 's', label = 'AUTOUGH2')
    plt.xlabel(field_name + ' (' + unit + ')')
    plt.ylabel('z (m)')
    plt.title(field_name + ' in column 0')
    img_filename_base = '_'.join((model_name, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(slice_test.mSuite.runs[run_index].basePath,
                                slice_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.legend(loc = 'best')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    slice_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], slice_test.outputPathBase):
    report_filename = os.path.join(slice_test.outputPathBase,
                     "%s-report.%s" % (slice_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(slice_test, mResults, rGen, report_filename)
    html_filename = os.path.join(slice_test.outputPathBase,
                     "%s-report.%s" % (slice_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
