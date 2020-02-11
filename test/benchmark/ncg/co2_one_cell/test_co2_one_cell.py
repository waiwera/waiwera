"""
One-cell CO2 test from O'Sullivan et al. (1985)
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import HistoryDataResult
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from credo.systest import HistoryWithinTolTC

from mulgrids import mulgrid

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'

import numpy as np
from docutils.core import publish_file

parser = argparse.ArgumentParser()
parser.add_argument("-np", type = int, default = 1, help = "number of processes") # not used
parser.add_argument("-d", "--docker", action = "store_true",
                    help = "run via Docker (waiwera-dkr)")
args = parser.parse_args()
num_procs = 1
mpi = False
simulator = 'waiwera-dkr -np %d' % num_procs if args.docker else 'waiwera'

model_name = 'co2_one_cell'

AUTOUGH2_FIELDMAP = {
    'Pressure': 'Pressure',
    'Temperature': 'Temperature',
    'Vapour saturation': 'Gas saturatio',
    'Enthalpy': 'Enthalpy'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Enthalpy': 'source_enthalpy'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)

run_name = 'run'
run_index = 0
obs_cell_index = 0
source_index = 0

test_fields = ['Pressure', 'Temperature', 'Vapour saturation']
test_source_fields = ['Enthalpy']
field_scale = {'Pressure': 1.e5, 'Temperature': 1.,
               'Vapour saturation': 1., 'Enthalpy': 1.e3}
field_unit = {'Pressure': 'bar', 'Temperature': '$^{\circ}$C',
               'Vapour saturation': '', 'Enthalpy': 'kJ/kg'}

digitised_test_fields = {'Pressure'}
digitised_simulators = ['MULKOM']

co2_one_cell_test = SciBenchmarkTest(model_name + "_test", nproc = num_procs)
co2_one_cell_test.description = """One-cell CO2 test from O'Sullivan et al. (1985), with Corey relative permeability curves and 30 bar initial CO2 partial pressure (fig 5).
"""

run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                          fieldname_map = WAIWERA_FIELDMAP,
                          simulator = simulator,
                          basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = num_procs
co2_one_cell_test.mSuite.addRun(model_run, run_name)

co2_one_cell_test.setupEmptyTestCompsList()

run_base_name = model_name
run_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", run_filename,
                                fieldname_map = AUTOUGH2_FIELDMAP,
                                 geo_filename = t2geo_filename)

co2_one_cell_test.addTestComp(run_index, "AUTOUGH2 history",
                      HistoryWithinTolTC(fieldsToTest = test_fields,
                                         defFieldTol = 1.e-3,
                                         expected = AUTOUGH2_result,
                                         testCellIndex = obs_cell_index))

co2_one_cell_test.addTestComp(run_index, "AUTOUGH2 source history",
                      HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                         defFieldTol = 1.e-3,
                                         expected = AUTOUGH2_result,
                                         testSourceIndex = source_index))

digitised_result = {}
for sim in digitised_simulators:
    data = {}
    for field_name in digitised_test_fields:
        data_filename = '_'.join((model_name, field_name, sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        data[field_name, obs_cell_index] = np.loadtxt(data_filename)
    digitised_result[sim] = HistoryDataResult(sim, data)

for sim in digitised_simulators:
    co2_one_cell_test.addTestComp(run_index, ' '.join((sim, field_name)),
                              HistoryWithinTolTC(fieldsToTest = digitised_test_fields,
                                                 defFieldTol = 2.e-2,
                                                 expected = digitised_result[sim],
                                                 testCellIndex = obs_cell_index,
                                                 orthogonalError = True))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = co2_one_cell_test.runTest(jrunner, createReports = True)

symbol = {"MULKOM": 'o'}

for field_name in digitised_test_fields:

    t, var = co2_one_cell_test.mSuite.resultsList[run_index].\
             getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.plot(t, var / field_scale[field_name], '-', label = 'Waiwera', zorder = 3)

    t, var = AUTOUGH2_result.getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.plot(t, var / field_scale[field_name], 's', label = 'AUTOUGH2', zorder = 2)

    for sim in digitised_simulators:
        result = digitised_result[sim]
        t, var = result.getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t, var / field_scale[field_name], symbol[sim], label = sim)

    plt.xlabel('time (s)')
    plt.ylabel(field_name + ' (' + field_unit[field_name] + ')')
    plt.legend(loc = 'best')
    img_filename_base = '_'.join((model_name, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(co2_one_cell_test.mSuite.runs[run_index].basePath,
                                co2_one_cell_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    co2_one_cell_test.mSuite.analysisImages.append(img_filename + '.png')

for field_name in test_source_fields:

    t, var = co2_one_cell_test.mSuite.resultsList[run_index].\
             getFieldHistoryAtSource(field_name, source_index)
    plt.plot(t, var / field_scale[field_name], '-', label = 'Waiwera', zorder = 3)

    t, var = AUTOUGH2_result.getFieldHistoryAtSource(field_name, source_index)
    plt.plot(t, var / field_scale[field_name], 's', label = 'AUTOUGH2', zorder = 2)

    plt.xlabel('time (s)')
    plt.ylabel(field_name + ' (' + field_unit[field_name] + ')')
    plt.legend(loc = 'best')
    img_filename_base = '_'.join((model_name, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(co2_one_cell_test.mSuite.runs[run_index].basePath,
                                co2_one_cell_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    co2_one_cell_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], co2_one_cell_test.outputPathBase):
    report_filename = os.path.join(co2_one_cell_test.outputPathBase,
                     "%s-report.%s" % (co2_one_cell_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(co2_one_cell_test, mResults, rGen, report_filename)
    html_filename = os.path.join(co2_one_cell_test.outputPathBase,
                     "%s-report.%s" % (co2_one_cell_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
