"""
Tests for sources on deliverability
"""

import os

from credo.systest import SciBenchmarkTest
from credo.systest import FieldWithinTolTC
from credo.systest import HistoryWithinTolTC

from credo.jobrunner import SimpleJobRunner
from credo.t2model import T2ModelRun, T2ModelResult
from credo.supermodel import SuperModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from mulgrids import mulgrid
import matplotlib.pyplot as plt
import numpy as np

from docutils.core import publish_file

model_name = 'deliv'

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
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

num_procs = 1

run_names = ['delv', 'delg_flow', 'delg_pi_table', 'delg_pwb_table', 'delg_limit', 'delt', 'delw']
test_fields = ["Pressure", "Temperature", "Vapour saturation"]
plot_fields = ["Pressure", "Vapour saturation"]

geo = mulgrid(t2geo_filename)
map_out_atm = range(geo.num_atmosphere_blocks, geo.num_blocks)

deliverability_test = SciBenchmarkTest("Deliverability", nproc = num_procs)
deliverability_test.description = """Tests sources on deliverability"""

for run_index, run_name in enumerate(run_names):

    run_base_name = '_'.join((model_name, run_name))
    run_filename = run_base_name + '.json'
    model_run = SuperModelRun(run_name, run_filename,
                              fieldname_map = SUPER_FIELDMAP,
                              simulator = 'supermodel',
                              basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = num_procs
    deliverability_test.mSuite.addRun(model_run, run_name)

deliverability_test.setupEmptyTestCompsList()

for run_index, run_name in enumerate(run_names):

    run_base_name = '_'.join((model_name, run_name))
    run_filename = os.path.join(model_dir, run_base_name + ".listing")
    reference_result = T2ModelResult("aut2", run_filename,
                                     geo_filename = t2geo_filename,
                                     ordering_map = map_out_atm,
                                     fieldname_map = AUT2_FIELDMAP)
    deliverability_test.addTestComp(run_index, "final errors",
                                        FieldWithinTolTC(fieldsToTest = test_fields,
                                                         defFieldTol = 5.e-3,
                                                         expected = reference_result,
                                                         testOutputIndex = -1))
    deliverability_test.addTestComp(run_index, "time history",
                                    HistoryWithinTolTC(fieldsToTest = test_fields,
                                                       defFieldTol = 1.e-2,
                                                       expected = reference_result,
                                                       testCellIndex = 4))
    
jrunner = SimpleJobRunner(mpi = True)
testResult, mResults = deliverability_test.runTest(jrunner, createReports = True)

# plots:
elevations = [layer.centre for layer in geo.layerlist[1:]]

for run_index, run_name in enumerate(run_names):

    run_base_name = '_'.join((model_name, run_name))
    tc_name = "time history"
    t = deliverability_test.testComps[run_index][tc_name].times
    for field_name in plot_fields:
        var = np.array(deliverability_test.testComps[run_index][tc_name].fieldErrors[field_name])
        plt.semilogx(t, var, '-o')
        plt.xlabel('time (s)')
        plt.ylabel(field_name + ' error')
        plt.title(run_name + ' ' + tc_name)
        img_filename_base = '_'.join((run_base_name, tc_name, 'error', field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(deliverability_test.mSuite.runs[run_index].basePath,
                                    deliverability_test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.savefig(img_filename)
        plt.clf()
        deliverability_test.mSuite.analysisImages.append(img_filename)

    tc_name = "final errors"
    for field_name in plot_fields:
        var = np.array(deliverability_test.testComps[run_index][tc_name].fieldErrors[field_name])
        plt.plot(var, elevations, 'o-')
        plt.ylabel('elevation (m)')
        plt.xlabel(field_name + ' error')
        plt.title(run_name + ' ' + tc_name)
        img_filename_base = '_'.join((run_base_name, tc_name, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(deliverability_test.mSuite.runs[run_index].basePath,
                                    deliverability_test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.savefig(img_filename)
        plt.clf()
        deliverability_test.mSuite.analysisImages.append(img_filename)

# generate report:

for rGen in getGenerators(["RST"], deliverability_test.outputPathBase):
    report_filename = os.path.join(deliverability_test.outputPathBase,
                     "%s-report.%s" % (deliverability_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(deliverability_test, mResults, rGen, report_filename)
    html_filename = os.path.join(deliverability_test.outputPathBase,
                     "%s-report.%s" % (deliverability_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")

