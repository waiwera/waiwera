"""
Model Intercomparison Study problem 6
"""

import os

from credo.systest import SciBenchmarkTest
from credo.systest import FieldWithinTolTC
from credo.systest import HistoryWithinTolTC

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import ModelResult
from credo.t2model import T2ModelRun, T2ModelResult
from credo.supermodel import SuperModelRun

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

model_name = 'problem6'

AUT2_FIELDMAP = {
    'Pressure': 'Pressure',
    'Temperature': 'Temperature',
    'Vapour saturation': 'Vapour saturation',
}
SUPER_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

num_procs = 1

run_name = 'run'
run_index = 0
test_fields = ["Pressure", "Temperature", "Vapour saturation"]
plot_fields = test_fields
digitised_test_fields = ["Pressure", "Vapour saturation"]
digitised_simulators = ["LBL", "S-Cubed"]    

geo = mulgrid(t2geo_filename)
map_out_atm = range(geo.num_atmosphere_blocks, geo.num_blocks)

problem6_test = SciBenchmarkTest(model_name + "_test", nproc = num_procs)
problem6_test.description = """Model Intercomparison Study problem 6"""

obspt = 'production'
obs_position = np.array([500., 400., -1050.])
obs_blk = geo.block_name_containing_point(obs_position)
obs_cell_index = geo.block_name_index[obs_blk] - geo.num_atmosphere_blocks

run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = SuperModelRun(run_name, run_filename,
                          fieldname_map = SUPER_FIELDMAP,
                          simulator = 'waiwera',
                          basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = num_procs
problem6_test.mSuite.addRun(model_run, run_name)

problem6_test.setupEmptyTestCompsList()
digitised_result = {}

run_base_name = model_name
run_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", run_filename,
                                 geo_filename = t2geo_filename,
                                 ordering_map = map_out_atm,
                                 fieldname_map = AUT2_FIELDMAP)

problem6_test.addTestComp(run_index, "AUTOUGH2 " + obspt + " well",
                      HistoryWithinTolTC(fieldsToTest = test_fields,
                                         defFieldTol = 1.e-3,
                                         expected = AUTOUGH2_result,
                                         testCellIndex = obs_cell_index))

for field_name in digitised_test_fields:
    for sim in digitised_simulators:
        data_filename = '_'.join((model_name, obspt, field_name, sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        result = DigitisedHistoryResult(sim, data_filename,
                                               field_name, obs_cell_index)

        digitised_result[obspt, field_name, sim] = result
        problem6_test.addTestComp(run_index, ' '.join((sim, field_name, obspt + ' well')),
                                  HistoryWithinTolTC(fieldsToTest = [field_name],
                                                     defFieldTol = 1.e-2,
                                                     expected = result,
                                                     testCellIndex = obs_cell_index,
                                                     orthogonalError = True))

jrunner = SimpleJobRunner(mpi = True)
testResult, mResults = problem6_test.runTest(jrunner, createReports = True)

# plots:
elevations = [layer.centre for layer in geo.layerlist[1:]]
scale = {"Pressure": 1.e5, "Temperature": 1., "Vapour saturation": 1.}
unit = {"Pressure": "bar", "Temperature": "$^{\circ}$C", "Vapour saturation": ""}
symbol = {"LBL": 's', "S-Cubed": 'o'}
yr = 365. * 24. * 60. * 60.

# plot time history results at well:
tc_name = "AUTOUGH2 " + obspt + " well"

for field_name in digitised_test_fields:
    t = problem6_test.testComps[run_index][tc_name].times
    var = problem6_test.mSuite.resultsList[run_index].getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.plot(t / yr, var / scale[field_name], '-', label = 'Waiwera')
    for sim in digitised_simulators:
        result = digitised_result[obspt, field_name, sim]
        t = result.getTimes()
        var = result.getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t / yr, var / scale[field_name], symbol[sim], label = sim)
    plt.xlabel('time (years)')
    plt.ylabel(field_name + ' (' + unit[field_name] + ')')
    plt.legend()
    plt.title(' '.join((model_name, field_name.lower(),
                        'results at', obspt, 'well')))
    img_filename_base = '_'.join((model_name, obspt, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem6_test.mSuite.runs[run_index].basePath,
                                problem6_test.mSuite.outputPathBase,
                                img_filename_base + '.png')
    plt.savefig(img_filename)
    plt.clf()
    problem6_test.mSuite.analysisImages.append(img_filename)

t = problem6_test.testComps[run_index][tc_name].times
for field_name in plot_fields:
    var = np.array(problem6_test.testComps[run_index][tc_name].fieldErrors[field_name])
    plt.plot(t / yr, var, '-o')
    plt.xlabel('time (years)')
    plt.ylabel(field_name + ' error')
    plt.title(' '.join((model_name, 'comparison with AUTOUGH2 at',
                       obspt, 'well')))
    img_filename_base = '_'.join((model_name, tc_name, 'error', field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem6_test.mSuite.runs[run_index].basePath,
                                problem6_test.mSuite.outputPathBase,
                                img_filename_base + '.png')
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

