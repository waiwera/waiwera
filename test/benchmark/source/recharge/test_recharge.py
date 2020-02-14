"""
Tests for recharge sources
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

model_name = 'recharge'

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Generation rate': 'source_rate',
    'Enthalpy': 'source_enthalpy'
}

model_dir = './run'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

run_names = ['outflow']
test_fields = ["Pressure", "Temperature"]
test_source_fields = ["Generation rate", "Enthalpy"]
plot_fields = ["Pressure"]
plot_source_fields = ["Generation rate"]
scale = {"Generation rate": 1.}
unit = {"Generation rate": "kg/s"}

geo = mulgrid(t2geo_filename)
map_out_atm = list(range(geo.num_atmosphere_blocks, geo.num_blocks))

test = SciBenchmarkTest("recharge_test", nproc = args.np)
test.description = """Tests recharge sources"""

source_index = 0

for run_index, run_name in enumerate(run_names):

    run_base_name = '_'.join((model_name, run_name))
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                              fieldname_map = WAIWERA_FIELDMAP,
                              simulator = simulator,
                              basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = args.np
    test.mSuite.addRun(model_run, run_name)

test.setupEmptyTestCompsList()

for run_index, run_name in enumerate(run_names):

    run_base_name = '_'.join((model_name, run_name))
    run_filename = os.path.join(model_dir, run_base_name + ".listing")
    reference_result = T2ModelResult("aut2", run_filename,
                                     geo_filename = t2geo_filename,
                                     ordering_map = map_out_atm)
    test.addTestComp(run_index, "final errors",
                                        FieldWithinTolTC(fieldsToTest = test_fields,
                                                         defFieldTol = 1.e-4,
                                                         expected = reference_result,
                                                         testOutputIndex = -1))
    test.addTestComp(run_index, "time history LH end",
                                    HistoryWithinTolTC(fieldsToTest = test_fields,
                                                       defFieldTol = 1.e-3,
                                                       expected = reference_result,
                                                       testCellIndex = 0))
    test.addTestComp(run_index, "time history RH end",
                                    HistoryWithinTolTC(fieldsToTest = test_fields,
                                                       defFieldTol = 1.e-3,
                                                       expected = reference_result,
                                                       testCellIndex = -1))
    test.addTestComp(run_index, "source",
                     HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                        defFieldTol = 1.e-2,
                                        expected = reference_result,
                                        testSourceIndex = source_index))
    
jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = test.runTest(jrunner, createReports = True)

# plots:
x = [col.centre[0] for col in geo.columnlist]

for run_index, run_name in enumerate(run_names):

    run_base_name = '_'.join((model_name, run_name))

    tc_name = "source"
    for field_name in plot_source_fields:

        t, var = test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtSource(field_name, source_index)
        plt.semilogx(t[1:], var[1:] / scale[field_name], '-', label = 'Waiwera')

        t, var = reference_result.getFieldHistoryAtSource(field_name, source_index)
        plt.semilogx(t[1:], var[1:] / scale[field_name], '+', label = 'AUTOUGH2')

        plt.xlabel('time (s)')
        plt.ylabel(field_name + ' (' + unit[field_name] + ')')
        plt.legend(loc = 'lower right')
        plt.title(' '.join((run_name, field_name.lower())))
        img_filename_base = '_'.join((model_name, run_name, tc_name, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(test.mSuite.runs[run_index].basePath,
                                    test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename)
        plt.clf()
        test.mSuite.analysisImages.append(img_filename)

    # plot time history results at end cells:
    for field_name in plot_fields:
        for cell in [1, -1]:
            t, var = test.mSuite.resultsList[run_index].\
                     getFieldHistoryAtCell(field_name, cell)
            plt.semilogx(t, var, '-o', label = 'Waiwera cell ' + str(cell))
            t, var = reference_result.getFieldHistoryAtCell(field_name, cell)
            plt.semilogx(t, var, 's', label = 'AUTOUGH2 cell ' + str(cell))
        plt.xlabel('t (s)')
        plt.ylabel(field_name)
        plt.legend()
        plt.title(run_name + ' results')
        img_filename_base = '_'.join((run_base_name, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(test.mSuite.runs[run_index].basePath,
                                    test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename)
        plt.clf()
        test.mSuite.analysisImages.append(img_filename)

# generate report:

for rGen in getGenerators(["RST"], test.outputPathBase):
    report_filename = os.path.join(test.outputPathBase,
                     "%s-report.%s" % (test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(test, mResults, rGen, report_filename)
    html_filename = os.path.join(test.outputPathBase,
                     "%s-report.%s" % (test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")

