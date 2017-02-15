"""
Model Intercomparison Study problem 1
"""

import os

from credo.systest import SciBenchmarkTest
from credo.systest import FieldWithinTolTC
from credo.systest import HistoryWithinTolTC

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import ModelResult
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from mulgrids import mulgrid
import matplotlib.pyplot as plt
import numpy as np

from docutils.core import publish_file

class DigitisedHistoryResult(ModelResult):
    """Digitised results for a field at a single cell."""
    def __init__(self, modelName, fileName, field, cellIndex, ordering_map=None,
                 fieldname_map=None):
        from os.path import dirname
        super(DigitisedHistoryResult, self).__init__(modelName, dirname(fileName),
                                            ordering_map=ordering_map,
                                            fieldname_map=fieldname_map)
        self.field = field
        self.cellIndex = cellIndex
        self.data = np.loadtxt(fileName)
    def _getTimes(self): return self.data[:, 0]
    def _getFieldHistoryAtCell(self, field, cellIndex):
        if field == self.field and cellIndex == self.cellIndex:
            return self.data[:, 1]
        else: return None

model_name = 'problem1'

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature'
}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

num_procs = 1

run_name = 'run'
run_index = 0
test_fields = ["Pressure", "Temperature"]
plot_fields = test_fields[:2]
digitised_test_fields = ["Temperature"]
digitised_simulators = ["S-Cubed", "GeoTrans"]

geo = mulgrid(t2geo_filename)
map_out_atm = range(geo.num_atmosphere_blocks, geo.num_blocks)

problem1_test = SciBenchmarkTest(model_name + "_test", nproc = num_procs)
problem1_test.description = """Model Intercomparison Study problem 1
(radial Avdonin problem)"""

obspt = 'r = 37.5 m'
obs_position = np.array([37.5, 0., -50.])
obs_blk = geo.block_name_containing_point(obs_position)
obs_cell_index = geo.block_name_index[obs_blk] - geo.num_atmosphere_blocks

run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                          fieldname_map = WAIWERA_FIELDMAP,
                          simulator = 'waiwera',
                          basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = num_procs
problem1_test.mSuite.addRun(model_run, run_name)

problem1_test.setupEmptyTestCompsList()
digitised_result = {}

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

for field_name in digitised_test_fields:
    for sim in digitised_simulators:
        data_filename = '_'.join((model_name, field_name, sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        result = DigitisedHistoryResult(sim, data_filename,
                                               field_name, obs_cell_index)

        digitised_result[obspt, field_name, sim] = result
        problem1_test.addTestComp(run_index, ' '.join((sim, field_name, obspt)),
                                  HistoryWithinTolTC(fieldsToTest = [field_name],
                                                     defFieldTol = 2.e-2,
                                                     expected = result,
                                                     testCellIndex = obs_cell_index,
                                                     orthogonalError = True,
                                                     logx = True))

jrunner = SimpleJobRunner(mpi = True)
testResult, mResults = problem1_test.runTest(jrunner, createReports = True)

# plots:
scale = {"Pressure": 1.e5, "Temperature": 1.}
unit = {"Pressure": "bar", "Temperature": "$^{\circ}$C"}
symbol = {"GeoTrans": 's', "S-Cubed": 'o'}

# plot time history results at r = 37.5 m:
tc_name = "AUTOUGH2 history at " + obspt

for field_name in digitised_test_fields:

    t = problem1_test.testComps[run_index][tc_name].times
    var = problem1_test.mSuite.resultsList[run_index].getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.semilogx(t, var / scale[field_name], '-', label = 'Waiwera')

    t = AUTOUGH2_result.getTimes()
    var = AUTOUGH2_result.getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.semilogx(t, var / scale[field_name], '+', label = 'AUTOUGH2')

    for sim in digitised_simulators:
        result = digitised_result[obspt, field_name, sim]
        t = result.getTimes()
        var = result.getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.semilogx(t, var / scale[field_name], symbol[sim], label = sim)
    plt.xlabel('time (s)')
    plt.ylabel(field_name + ' (' + unit[field_name] + ')')
    plt.legend(loc = 'best')
    plt.title(' '.join((model_name, field_name.lower(),
                        'results at', obspt)))
    img_filename_base = '_'.join((model_name, obspt, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem1_test.mSuite.runs[run_index].basePath,
                                problem1_test.mSuite.outputPathBase,
                                img_filename_base + '.png')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename)
    plt.clf()
    problem1_test.mSuite.analysisImages.append(img_filename)

t = problem1_test.testComps[run_index][tc_name].times
for field_name in plot_fields:
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
                                img_filename_base + '.png')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename)
    plt.clf()
    problem1_test.mSuite.analysisImages.append(img_filename)

# plot temperature profile at end time:
tc_name = "AUTOUGH2 at t = 1.e9 s"
r = np.array([col.centre[0] for col in geo.columnlist])
for field_name in plot_fields:
    var = problem1_test.mSuite.resultsList[run_index].getFieldAtOutputIndex(field_name, -1) /scale[field_name]
    plt.plot(r, var, '-', label = 'Waiwera')
    var = AUTOUGH2_result.getFieldAtOutputIndex(field_name, -1) / scale[field_name]
    plt.plot(r, var, '+', label = 'AUTOUGH2')
    plt.xlabel('radius (m)')
    plt.ylabel(field_name + ' (' + unit[field_name] + ')')
    plt.title(' '.join((model_name, 'comparison with', tc_name)))
    img_filename_base = '_'.join((model_name, tc_name, 'comparison', field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem1_test.mSuite.runs[run_index].basePath,
                                problem1_test.mSuite.outputPathBase,
                                img_filename_base + '.png')
    plt.legend(loc = 'best')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename)
    plt.clf()
    problem1_test.mSuite.analysisImages.append(img_filename)

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

