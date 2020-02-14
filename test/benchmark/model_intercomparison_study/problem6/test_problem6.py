"""
Model Intercomparison Study problem 6
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest
from credo.systest import FieldWithinTolTC
from credo.systest import HistoryWithinTolTC

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import ModelResult, HistoryDataResult
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

model_name = 'problem6'

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Enthalpy': 'source_enthalpy'
}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

run_name = 'run'
run_index = 0
test_fields = ["Pressure", "Temperature", "Vapour saturation"]
plot_fields = test_fields
digitised_test_fields = ["Pressure", "Vapour saturation"]
test_source_fields = ["Enthalpy"]
digitised_simulators = ["LBL", "S-Cubed"]    

geo = mulgrid(t2geo_filename)
map_out_atm = list(range(geo.num_atmosphere_blocks, geo.num_blocks))

problem6_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
problem6_test.description = """Model Intercomparison Study problem 6"""

obspt = 'production'
obs_position = np.array([500., 400., -1050.])
obs_blk = geo.block_name_containing_point(obs_position)
obs_cell_index = geo.block_name_index[obs_blk] - geo.num_atmosphere_blocks
source_index = 0

run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                          fieldname_map = WAIWERA_FIELDMAP,
                          simulator = simulator,
                          basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = args.np
problem6_test.mSuite.addRun(model_run, run_name)

problem6_test.setupEmptyTestCompsList()

run_base_name = model_name
run_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", run_filename,
                                 geo_filename = t2geo_filename,
                                 ordering_map = map_out_atm)

problem6_test.addTestComp(run_index, "AUTOUGH2 " + obspt + ' well',
                          HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                             defFieldTol = 2.e-2,
                                             expected = AUTOUGH2_result,
                                             testSourceIndex = source_index))
digitised_result = {}
digitised_source_result = {}
source_index = 0
for sim in digitised_simulators:
    data = {}
    for field_name in digitised_test_fields:
        data_filename = '_'.join((model_name, obspt, field_name, sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        data[field_name, obs_cell_index] = np.loadtxt(data_filename)
    digitised_result[sim] = HistoryDataResult(sim, data)
    data = {}
    for field_name in test_source_fields:
        data_filename = '_'.join((model_name, obspt, field_name, sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        data[field_name, source_index] = np.loadtxt(data_filename)
    digitised_source_result[sim] = HistoryDataResult(sim, data)

sim_tol = {"LBL": 0.075, "S-Cubed": 0.015}

for sim in digitised_simulators:
    problem6_test.addTestComp(run_index, ' '.join((sim, obspt + ' well')),
                              HistoryWithinTolTC(fieldsToTest = digitised_test_fields,
                                                 defFieldTol = sim_tol[sim],
                                                 expected = digitised_result[sim],
                                                 testCellIndex = obs_cell_index,
                                                 orthogonalError = True))
    problem6_test.addTestComp(run_index, ' '.join((sim, obspt + ' well discharge')),
                              HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                                 defFieldTol = sim_tol[sim],
                                                 expected = digitised_source_result[sim],
                                                 testSourceIndex = source_index,
                                                 orthogonalError = True))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = problem6_test.runTest(jrunner, createReports = True)

# plots:
elevations = [layer.centre for layer in geo.layerlist[1:]]
scale = {"Pressure": 1.e5, "Temperature": 1., "Vapour saturation": 1.,
         "Enthalpy": 1.e3}
unit = {"Pressure": "bar", "Temperature": "$^{\circ}$C", "Vapour saturation": "",
        "Enthalpy": "kJ/kg"}
symbol = {"LBL": 's', "S-Cubed": 'o'}
yr = 365.25 * 24. * 60. * 60.

# plot time history results at well:

for field_name in digitised_test_fields:

    t, var = problem6_test.mSuite.resultsList[run_index].\
             getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.plot(t / yr, var / scale[field_name], '-', label = 'Waiwera')

    t, var = AUTOUGH2_result.getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.plot(t / yr, var / scale[field_name], '+', label = 'AUTOUGH2')

    for sim in digitised_simulators:
        result = digitised_result[sim]
        t, var = result.getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t / yr, var / scale[field_name], symbol[sim], label = sim)
    plt.xlabel('time (years)')
    plt.ylabel(field_name + ' (' + unit[field_name] + ')')
    plt.legend(loc = 'best')
    plt.title(' '.join((model_name, field_name.lower(),
                        'results at', obspt, 'well')))
    img_filename_base = '_'.join((model_name, obspt, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem6_test.mSuite.runs[run_index].basePath,
                                problem6_test.mSuite.outputPathBase,
                                img_filename_base + '.png')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename)
    plt.clf()
    problem6_test.mSuite.analysisImages.append(img_filename)

for field_name in test_source_fields:

    t, var = problem6_test.mSuite.resultsList[run_index].\
             getFieldHistoryAtSource(field_name, source_index)
    plt.plot(t / yr, var / scale[field_name], '-', label = 'Waiwera')

    t, var = AUTOUGH2_result.getFieldHistoryAtSource(field_name, source_index)
    plt.plot(t / yr, var / scale[field_name], '+', label = 'AUTOUGH2')

    for sim in digitised_simulators:
        result = digitised_source_result[sim]
        t, var = result.getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t / yr, var / scale[field_name], symbol[sim], label = sim)
    plt.xlabel('time (years)')
    plt.ylabel(field_name + ' (' + unit[field_name] + ')')
    plt.legend(loc = 'best')
    plt.title(' '.join((model_name, field_name.lower(),
                        'results at', obspt, 'well')))
    img_filename_base = '_'.join((model_name, obspt, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem6_test.mSuite.runs[run_index].basePath,
                                problem6_test.mSuite.outputPathBase,
                                img_filename_base + '.png')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename)
    plt.clf()
    problem6_test.mSuite.analysisImages.append(img_filename)

# generate report:

for rGen in getGenerators(["RST"], problem6_test.outputPathBase):
    report_filename = os.path.join(problem6_test.outputPathBase,
                     "%s-report.%s" % (problem6_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(problem6_test, mResults, rGen, report_filename)
    html_filename = os.path.join(problem6_test.outputPathBase,
                     "%s-report.%s" % (problem6_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")

