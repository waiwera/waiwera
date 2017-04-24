"""
TOUGH2 infiltration test from TOUGH user guide
"""

import os

from credo.systest import SciBenchmarkTest

from credo.jobrunner import SimpleJobRunner
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun
from credo.modelresult import DigitisedOneDFieldResult

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from credo.systest import FieldWithinTolTC, OneDSolutionWithinTolTC

from mulgrids import mulgrid

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'

import numpy as np
from docutils.core import publish_file
    
model_name = 'infiltration'

def liquid_sat(mResult, index):
    sv = mResult.getFieldAtOutputIndex('Gas saturati', index)
    return 1. - sv
    
AUTOUGH2_FIELDMAP = {
    'Liquid saturation': liquid_sat}

WAIWERA_FIELDMAP = {
    'Liquid saturation': 'fluid_liquid_saturation'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
map_out_bdy = range(0, geo.num_blocks)

num_procs = 1
run_index = 0
run_name = 'run'
run_base_name = model_name
output_indices = [1, 2, 3]

test_fields = ['Liquid saturation']
digitised_test_fields = test_fields
digitised_simulators = ["semi-analytical"]

infiltration_test = SciBenchmarkTest(model_name + "_test", nproc = num_procs)
infiltration_test.description = """1-D horizontal infiltration test from TOUGH User Guide
"""

run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                            fieldname_map = WAIWERA_FIELDMAP,
                            simulator = 'waiwera',
                            basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = num_procs
infiltration_test.mSuite.addRun(model_run, run_name)

infiltration_test.setupEmptyTestCompsList()

results_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", results_filename,
                                geo_filename = t2geo_filename,
                                fieldname_map = AUTOUGH2_FIELDMAP,
                                ordering_map = map_out_bdy)
for output_index in output_indices:
    infiltration_test.addTestComp(run_index, "AUTOUGH2",
                                  FieldWithinTolTC(fieldsToTest = test_fields,
                                                   defFieldTol = 1.e-4,
                                                   expected = AUTOUGH2_result,
                                                   testOutputIndex = output_index))

digitised_result = {}
digitised_times = [864, 5184, 9504]
digitised_output_index = [1, 2, 3]
xmax_all = 0.
for sim in digitised_simulators:
    for field_name in digitised_test_fields:
        for time, output_index in zip(digitised_times, digitised_output_index):
            data_filename = '_'.join([sim, str(time)])
            data_filename = data_filename.lower().replace(' ', '_')
            data_filename = os.path.join(data_dir, data_filename + '.dat')
            result = DigitisedOneDFieldResult(sim, data_filename, field_name, output_index)
            digitised_result[sim, output_index] = result
            xmax = np.max(result.getCoordinates())
            xmax_all = max(xmax_all, xmax)
            infiltration_test.addTestComp(run_index, ' '.join((sim, str(output_index))),
                                      OneDSolutionWithinTolTC(
                                          fieldsToTest = [field_name],
                                          defFieldTol = 5.e-2,
                                          expected = result,
                                          maxCoordinate = xmax,
                                          testOutputIndex = output_index))

jrunner = SimpleJobRunner(mpi = True)
testResult, mResults = infiltration_test.runTest(jrunner, createReports = True)

result = infiltration_test.mSuite.resultsList[run_index]
t = result.getTimes()
x = np.array([col.centre[0] for col in geo.columnlist])
ix = np.where(x <= xmax_all)

for field_name in test_fields:
    for output_index in output_indices:
        var = result.getFieldAtOutputIndex(field_name, output_index)
        tstr = ' t = %4.0f s' % t[output_index]
        plt.plot(x[ix], var[ix], 'o', label = 'Waiwera' + tstr)
        var = AUTOUGH2_result.getFieldAtOutputIndex(field_name, output_index)
        plt.plot(x[ix], var[ix], '+', label = 'AUTOUGH2' + tstr)
        for sim in digitised_simulators:
            dig_result = digitised_result[sim, output_index]
            xd = dig_result.getCoordinates()
            var = dig_result.getFieldAtOutputIndex(field_name, output_index)
            plt.plot(xd, var, '-', label = sim + tstr)
    plt.xlabel('x (m)')
    plt.ylabel(field_name)
    img_filename_base = '_'.join((model_name, run_name, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(infiltration_test.mSuite.runs[run_index].basePath,
                                infiltration_test.mSuite.outputPathBase,
                                img_filename_base + '.png')
    plt.legend(loc = 'best')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename)
    plt.clf()
    infiltration_test.mSuite.analysisImages.append(img_filename)

# generate report:

for rGen in getGenerators(["RST"], infiltration_test.outputPathBase):
    report_filename = os.path.join(infiltration_test.outputPathBase,
                     "%s-report.%s" % (infiltration_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(infiltration_test, mResults, rGen, report_filename)
    html_filename = os.path.join(infiltration_test.outputPathBase,
                     "%s-report.%s" % (infiltration_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
