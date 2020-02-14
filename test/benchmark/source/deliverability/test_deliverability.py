"""
Tests for sources on deliverability
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

model_name = 'deliv'

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Generation rate': 'source_rate',
    'Enthalpy': 'source_enthalpy'
}

model_dir = './run'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

run_names = ['delv', 'delg_flow', 'delg_pi_table', 'delg_pwb_table', 'delg_limit', 'delt', 'delw']
test_fields = ["Pressure", "Temperature", "Vapour saturation"]
test_source_fields = ["Generation rate", "Enthalpy"]
plot_fields = ["Pressure", "Vapour saturation"]
plot_source_fields = ["Generation rate"]
scale = {"Generation rate": 1.}
unit = {"Generation rate": "kg/s"}

geo = mulgrid(t2geo_filename)
map_out_atm = list(range(geo.num_atmosphere_blocks, geo.num_blocks))

deliverability_test = SciBenchmarkTest("deliverability_test", nproc = args.np)
deliverability_test.description = """Tests sources on deliverability"""

source_index = 0

for run_index, run_name in enumerate(run_names):

    run_base_name = '_'.join((model_name, run_name))
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                              fieldname_map = WAIWERA_FIELDMAP,
                              simulator = simulator,
                              basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = args.np
    deliverability_test.mSuite.addRun(model_run, run_name)

deliverability_test.setupEmptyTestCompsList()

reference_result = {}
for run_index, run_name in enumerate(run_names):

    run_base_name = '_'.join((model_name, run_name))
    run_filename = os.path.join(model_dir, run_base_name + ".listing")
    reference_result[run_name] = T2ModelResult("aut2", run_filename,
                                     geo_filename = t2geo_filename,
                                     ordering_map = map_out_atm)
    deliverability_test.addTestComp(run_index, "final errors",
                                        FieldWithinTolTC(fieldsToTest = test_fields,
                                                         defFieldTol = 5.e-3,
                                                         expected = reference_result[run_name],
                                                         testOutputIndex = -1))
    deliverability_test.addTestComp(run_index, "time history",
                                    HistoryWithinTolTC(fieldsToTest = test_fields,
                                                       defFieldTol = 1.e-2,
                                                       expected = reference_result[run_name],
                                                       testCellIndex = 4))
    deliverability_test.addTestComp(run_index, "source",
                                    HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                                       defFieldTol = 1.e-2,
                                                       expected = reference_result[run_name],
                                                       testSourceIndex = source_index))
    
jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = deliverability_test.runTest(jrunner, createReports = True)

# plots:
x = [col.centre[0] for col in geo.columnlist]

for run_index, run_name in enumerate(run_names):

    tc_name = "source"
    for field_name in plot_source_fields:

        t, var = deliverability_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t, var / scale[field_name], '-', label = 'Waiwera')

        t, var = reference_result[run_name].getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t, var / scale[field_name], '+', label = 'AUTOUGH2')

        plt.xlabel('time (s)')
        plt.ylabel(field_name + ' (' + unit[field_name] + ')')
        plt.legend(loc = 'best')
        plt.title(' '.join((run_name, field_name.lower())))
        img_filename_base = '_'.join((model_name, run_name, tc_name, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(deliverability_test.mSuite.runs[run_index].basePath,
                                    deliverability_test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.tight_layout(pad = 3.)
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

