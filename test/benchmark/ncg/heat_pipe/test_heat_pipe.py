"""
Radial heat pipe problem ('rhp') from TOUGH2 user guide
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import DigitisedOneDFieldResult
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from credo.systest import FieldWithinTolTC, OneDSolutionWithinTolTC

from mulgrids import mulgrid

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

model_name = 'heat_pipe'

AUTOUGH2_FIELDMAP = {
    'Pressure': 'Pressure',
    'Temperature': 'Temperature',
    'Vapour saturation': 'Gas saturati',
    'Vapour air mass fraction': 'Air gas mass'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Vapour air mass fraction': 'fluid_vapour_air_mass_fraction'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
min_radius, max_radius = 1.e-1, 400.

run_name = 'run'
run_index = 0

test_fields = {'Pressure', 'Temperature', 'Vapour saturation',
               'Vapour air mass fraction'}
field_scale = {'Pressure': 1.e5, 'Temperature': 1.,
               'Vapour saturation': 1., 'Vapour air mass fraction': 1.}
field_unit = {'Pressure': 'bar', 'Temperature': '$^{\circ}$C',
               'Vapour saturation': '', 'Vapour air mass fraction': ''}

digitised_test_fields = test_fields
digitised_simulators = ["TOUGH2"]

map_out_bdy = list(range(0, geo.num_blocks))

heat_pipe_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
heat_pipe_test.description = """Radial heat pipe problem,
from the TOUGH2 user guide. (problem 'rhp').
Here the TOUGH2 results are digitised from the user guide. There is some uncertainty about the
digitised TOUGH2 pressure results at small radii, which don't match the AUTOUGH2 results.
This may have resulted from the small pressure range compared with the plot scale."""

run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                          fieldname_map = WAIWERA_FIELDMAP,
                          simulator = simulator,
                          basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = args.np
heat_pipe_test.mSuite.addRun(model_run, run_name)

heat_pipe_test.setupEmptyTestCompsList()
digitised_result = {}

run_base_name = model_name
run_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", run_filename,
                                fieldname_map = AUTOUGH2_FIELDMAP,
                                 geo_filename = t2geo_filename,
                                 ordering_map = map_out_bdy)

heat_pipe_test.addTestComp(run_index, "AUTOUGH2",
                           FieldWithinTolTC(fieldsToTest = test_fields,
                                            defFieldTol = 5.e-3,
                                            expected = AUTOUGH2_result,
                                            testOutputIndex = -1))

t_final = AUTOUGH2_result.getTimes()[-1]
for sim in digitised_simulators:    
    for field_name in digitised_test_fields:
        data_filename = '_'.join((model_name, field_name, sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        result = DigitisedOneDFieldResult(sim, data_filename, field_name, -1)
        result.data[:,0] = np.exp(result.data[:,0]) * np.sqrt(t_final)
        result.data[:,1] *= field_scale[field_name]
        digitised_result[field_name, sim] = result
        heat_pipe_test.addTestComp(run_index, ' '.join((sim, field_name)),
                                  OneDSolutionWithinTolTC(
                                      fieldsToTest = [field_name],
                                      defFieldTol = 2.e-2,
                                      fieldTols = {'Pressure': 0.07},
                                      expected = result,
                                      maxCoordinate = max_radius,
                                      logCoordinate = True,
                                      testOutputIndex = -1))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = heat_pipe_test.runTest(jrunner, createReports = True)

# plots:
symbol = {"TOUGH2": 'o'}

tc_name = "AUTOUGH2"
outputIndex = -1

for field_name in test_fields:
    result = heat_pipe_test.mSuite.resultsList[run_index]
    r = np.array([col.centre[0] for col in geo.columnlist])
    var = result.getFieldAtOutputIndex(field_name, outputIndex) / field_scale[field_name]
    plt.semilogx(r, var, '-', label = 'Waiwera', zorder = 3)
    var = AUTOUGH2_result.getFieldAtOutputIndex(field_name,
                                                outputIndex) / field_scale[field_name]
    plt.semilogx(r, var, 's', label = 'AUTOUGH2', zorder = 2)
    for sim in digitised_simulators:
        result = digitised_result[field_name, sim]
        r = result.getCoordinates()
        var = result.getFieldAtOutputIndex(field_name, outputIndex)
        plt.semilogx(r, var / field_scale[field_name], symbol[sim], label = sim)
    plt.xlabel('radius (m)')
    plt.ylabel(field_name + ' (' + field_unit[field_name] + ')')
    plt.title(field_name)
    img_filename_base = '_'.join((model_name, tc_name, 'comparison', field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(heat_pipe_test.mSuite.runs[run_index].basePath,
                                heat_pipe_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.xlim([min_radius, max_radius])
    plt.legend(loc = 'center right')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    heat_pipe_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], heat_pipe_test.outputPathBase):
    report_filename = os.path.join(heat_pipe_test.outputPathBase,
                     "%s-report.%s" % (heat_pipe_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(heat_pipe_test, mResults, rGen, report_filename)
    html_filename = os.path.join(heat_pipe_test.outputPathBase,
                     "%s-report.%s" % (heat_pipe_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
