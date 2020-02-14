"""
Model Intercomparison Study problem 1
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest
from credo.systest import FieldWithinTolTC, HistoryWithinTolTC, OneDSolutionWithinTolTC

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import DigitisedOneDFieldResult, HistoryDataResult
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from mulgrids import mulgrid
import matplotlib.pyplot as plt
import numpy as np

from docutils.core import publish_file

parser = argparse.ArgumentParser()
parser.add_argument("-np", type = int, default = 1, help = "number of processes")
parser.add_argument("-d", "--docker", action = "store_true",
                    help = "run via Docker (waiwera-dkr)")
args = parser.parse_args()
mpi = args.np > 1 and not args.docker
simulator = 'waiwera-dkr -np %d' % args.np if args.docker else 'waiwera'

model_name = 'problem1'

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature'
}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

run_name = 'run'
run_index = 0
test_fields = ["Pressure", "Temperature"]
digitised_test_fields = ["Temperature"]
digitised_simulators = ["S-Cubed"]

geo = mulgrid(t2geo_filename)
map_out_atm = list(range(geo.num_atmosphere_blocks, geo.num_blocks))

problem1_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
problem1_test.description = """Model Intercomparison Study problem 1
(radial Avdonin problem)"""

obspt = 'r = 37.5 m'
obs_position = np.array([37.5, 0., -50.])
obs_blk = geo.block_name_containing_point(obs_position)
obs_cell_index = geo.block_name_index[obs_blk] - geo.num_atmosphere_blocks
max_radius = 500.

run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                          fieldname_map = WAIWERA_FIELDMAP,
                          simulator = simulator,
                          basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = args.np
problem1_test.mSuite.addRun(model_run, run_name)

problem1_test.setupEmptyTestCompsList()
digitised_time_result = {}
digitised_r_result = {}

run_base_name = model_name
run_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", run_filename,
                                 geo_filename = t2geo_filename,
                                 ordering_map = map_out_atm)


problem1_test.addTestComp(run_index, "AUTOUGH2 history at " + obspt,
                      HistoryWithinTolTC(fieldsToTest = test_fields,
                                         defFieldTol = 1.e-3,
                                         expected = AUTOUGH2_result,
                                         testCellIndex = obs_cell_index))

problem1_test.addTestComp(run_index, "AUTOUGH2 t = 1.e9 s",
        FieldWithinTolTC(fieldsToTest = ["Temperature"],
                         defFieldTol = 1.0e-4,
                         expected = AUTOUGH2_result,
                         testOutputIndex = -1))

digitised_result = {}
for sim in digitised_simulators:
    data = {}
    for field_name in digitised_test_fields:
        data_filename = '_'.join((model_name, field_name, 'time', sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        data[field_name, obs_cell_index] = np.loadtxt(data_filename)
    digitised_result[sim] = HistoryDataResult(sim, data)

for sim in digitised_simulators:
    problem1_test.addTestComp(run_index, ' '.join((sim, field_name, obspt)),
                              HistoryWithinTolTC(fieldsToTest = [field_name],
                                                 defFieldTol = 2.e-2,
                                                 expected = digitised_result[sim],
                                                 testCellIndex = obs_cell_index,
                                                 orthogonalError = True,
                                                 logx = True))
    
    for field_name in digitised_test_fields:
        data_filename = '_'.join((model_name, field_name, 'r', sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        result = DigitisedOneDFieldResult(sim, data_filename, field_name, -1)
        digitised_r_result[field_name, sim] = result
        problem1_test.addTestComp(run_index, ' '.join((sim, field_name)),
                                  OneDSolutionWithinTolTC(
                                      fieldsToTest = [field_name],
                                      defFieldTol = 2.e-2,
                                      expected = result,
                                      testOutputIndex = -1,
                                      maxCoordinate = max_radius,
                                      logCoordinate = True))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = problem1_test.runTest(jrunner, createReports = True)

# plots:
scale = {"Pressure": 1.e5, "Temperature": 1.}
unit = {"Pressure": "bar", "Temperature": "$^{\circ}$C"}
symbol = {"GeoTrans": 's', "S-Cubed": 'o'}

# plot time history results at r = 37.5 m:
tc_name = "AUTOUGH2 history at " + obspt

sim = 'analytical'
data = {}
for field_name in digitised_test_fields:
    data_filename = '_'.join((model_name, field_name, 'time', sim))
    data_filename = data_filename.lower().replace(' ', '_')
    data_filename = os.path.join(data_dir, data_filename + '.dat')
    data[field_name, obs_cell_index] = np.loadtxt(data_filename)
analytical_result = HistoryDataResult(sim, data)

for field_name in digitised_test_fields:

    t, var = problem1_test.mSuite.resultsList[run_index].\
             getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.semilogx(t, var / scale[field_name], '-', label = 'Waiwera', zorder = 4)

    t, var = AUTOUGH2_result.getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.semilogx(t, var / scale[field_name], 's', label = 'AUTOUGH2', zorder = 3)

    for sim in digitised_simulators:
        result = digitised_result[sim]
        t, var = result.getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.semilogx(t, var / scale[field_name], symbol[sim], label = sim)

    t, var = analytical_result.getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.semilogx(t, var / scale[field_name], 'k--', label = 'analytical', zorder = 1)
    plt.xlabel('time (s)')
    plt.ylabel(field_name + ' (' + unit[field_name] + ')')

    plt.legend(loc = 'best')
    plt.title(' '.join((model_name, field_name.lower(),
                        'results at', obspt)))
    img_filename_base = '_'.join((model_name, obspt, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem1_test.mSuite.runs[run_index].basePath,
                                problem1_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    problem1_test.mSuite.analysisImages.append(img_filename + '.png')

t = problem1_test.testComps[run_index][tc_name].times
for field_name in digitised_test_fields:
    var = np.array(problem1_test.testComps[run_index][tc_name].fieldErrors[field_name])
    plt.semilogx(t, var, '-o')
    plt.xlabel('time (s)')
    plt.ylabel(field_name + ' error')
    plt.title(' '.join((model_name, 'comparison with AUTOUGH2 at',
                       obspt)))
    img_filename_base = '_'.join((model_name, tc_name, 'error', field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem1_test.mSuite.runs[run_index].basePath,
                                problem1_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    problem1_test.mSuite.analysisImages.append(img_filename + '.png')

# plot temperature profile w.r.t. radius at end time:
tc_name = "AUTOUGH2 at t = 1.e9 s"
outputIndex = -1
r = np.array([col.centre[0] for col in geo.columnlist])
ir = np.where(r <= max_radius)
r = r[ir]

for field_name in digitised_test_fields:
    result = problem1_test.mSuite.resultsList[run_index]
    var = result.getFieldAtOutputIndex(field_name, outputIndex)[ir] /scale[field_name]
    plt.plot(r, var, '-', label = 'Waiwera', zorder = 4)
    var = AUTOUGH2_result.getFieldAtOutputIndex(field_name,
                                                outputIndex)[ir] / scale[field_name]
    plt.plot(r, var, 's', label = 'AUTOUGH2', zorder = 3)
    for sim in digitised_simulators:
        result = digitised_r_result[field_name, sim]
        r = result.getCoordinates()
        var = result.getFieldAtOutputIndex(field_name, outputIndex)
        plt.plot(r, var / scale[field_name], symbol[sim], label = sim)
    sim = 'analytical'
    data_filename = '_'.join((model_name, field_name, 'r', sim))
    data_filename = data_filename.lower().replace(' ', '_')
    data_filename = os.path.join(data_dir, data_filename + '.dat')
    result = DigitisedOneDFieldResult(sim, data_filename, field_name, -1)
    r = result.getCoordinates()
    var = result.getFieldAtOutputIndex(field_name, outputIndex)
    plt.plot(r, var / scale[field_name], 'k--', label = sim, zorder = 1)
    plt.xlabel('radius (m)')
    plt.ylabel(field_name + ' (' + unit[field_name] + ')')
    plt.title(' '.join((model_name, 'comparison with', tc_name)))
    img_filename_base = '_'.join((model_name, tc_name, 'comparison', field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem1_test.mSuite.runs[run_index].basePath,
                                problem1_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.legend(loc = 'upper left')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    problem1_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], problem1_test.outputPathBase):
    report_filename = os.path.join(problem1_test.outputPathBase,
                     "%s-report.%s" % (problem1_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(problem1_test, mResults, rGen, report_filename)
    html_filename = os.path.join(problem1_test.outputPathBase,
                     "%s-report.%s" % (problem1_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")

