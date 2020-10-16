"""
1-D single-phase liquid tracer test with Dirichlet upstream boundary condition
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest

from credo.jobrunner import SimpleJobRunner
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun
from credo.modelresult import DigitisedOneDFieldResult

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from credo.systest import FieldWithinTolTC
from credo.systest import HistoryWithinTolTC

from mulgrids import mulgrid

import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'

day = 24. * 60. * 60.
import numpy as np
from docutils.core import publish_file

parser = argparse.ArgumentParser()
parser.add_argument("-np", type = int, default = 1, help = "number of processes")
parser.add_argument("-d", "--docker", action = "store_true",
                    help = "run via Docker (waiwera-dkr)")
args = parser.parse_args()
mpi = args.np > 1 and not args.docker
simulator = 'waiwera-dkr -np %d' % args.np if args.docker else 'waiwera'

model_name = 'oned'
    
AUTOUGH2_FIELDMAP = {
    'Pressure': 'Pressure',
    'Tracer mass fraction': 'Tracer/liquid',
    'Tracer production rate': 'Tracer mass flow'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Tracer mass fraction': 'tracer_tracer',
    'Tracer production rate': 'source_tracer_flow'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
map_out_bdy = list(range(0, geo.num_blocks))
source_index = 0

run_names = ['single', 'two']

test_fields = ['Pressure', 'Tracer mass fraction']
test_source_fields = ['Tracer production rate']
field_unit = {'Pressure': 'bar', 'Tracer mass fraction': '-',
              'Tracer production rate': 'mg/s'}
field_scale = {'Pressure': 1.e5, 'Tracer mass fraction': 1,
               'Tracer production rate': 1.e-6}

tracer_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
tracer_test.description = """1-D single-phase and two-phase liquid tracer test cases
with Dirichlet upstream boundary condition
"""
for run_index, run_name in enumerate(run_names):
    run_base_name = model_name + '_' + run_name + '_phase'
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                                fieldname_map = WAIWERA_FIELDMAP,
                                simulator = simulator,
                                basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = args.np
    tracer_test.mSuite.addRun(model_run, run_name)

tracer_test.setupEmptyTestCompsList()
AUTOUGH2_result = {}

for run_index, run_name in enumerate(run_names):
    run_base_name = model_name + '_' + run_name + '_phase'
    results_filename = os.path.join(model_dir, run_base_name + ".listing")
    AUTOUGH2_result[run_name] = T2ModelResult("AUTOUGH2", results_filename,
                                              geo_filename = t2geo_filename,
                                              fieldname_map = AUTOUGH2_FIELDMAP,
                                              ordering_map = map_out_bdy)
    tracer_test.addTestComp(run_index, "AUTOUGH2",
                            FieldWithinTolTC(fieldsToTest = test_fields,
                                             defFieldTol = 1.e-3,
                                             absoluteErrorTol = 1.e-4,
                                             expected = AUTOUGH2_result[run_name],
                                             testOutputIndex = -1))
    tracer_test.addTestComp(run_index, 'AUTOUGH2 source',
                            HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                               defFieldTol = 1.e-3,
                                               expected = AUTOUGH2_result[run_name],
                                               testSourceIndex = source_index))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = tracer_test.runTest(jrunner, createReports = True)

x = np.array([col.centre[0] for col in geo.columnlist])

for run_index, run_name in enumerate(run_names):

    for field_name in test_fields:
        scale = field_scale[field_name]
        unit = field_unit[field_name]
        result = tracer_test.mSuite.resultsList[run_index]
        t = result.getTimes()
        var = result.getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(x, var, 'b-', label = 'Waiwera', zorder = 2)
        var = AUTOUGH2_result[run_name].getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(x, var, 'gs', label = 'AUTOUGH2', zorder = 1)
        plt.xlabel('x (m)')
        plt.ylabel('%s (%s)' % (field_name, unit))
        img_filename_base = '_'.join((model_name, run_name, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(tracer_test.mSuite.runs[run_index].basePath,
                                    tracer_test.mSuite.outputPathBase,
                                    img_filename_base)
        plt.legend(loc = 'best')
        plt.title('%s at time %2.0f days: %s phase' % (field_name, t[-1] / day, run_name))
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename + '.png', dpi = 300)
        plt.savefig(img_filename + '.pdf')
        plt.clf()
        tracer_test.mSuite.analysisImages.append(img_filename + '.png')

    for field_name in test_source_fields:
        t, var = tracer_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t / day, var / field_scale[field_name], 'b-', label = 'Waiwera')

        t, var = AUTOUGH2_result[run_name].getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t / day, var / field_scale[field_name], 'gs', label = 'AUTOUGH2')

        plt.xlabel('time (days)')
        plt.ylabel(field_name + ' (' + field_unit[field_name] + ')')
        plt.legend(loc = 'best')
        plt.title(field_name)
        plt.title('%s history: %s phase' % (field_name, run_name))
        img_filename_base = '_'.join((model_name, run_name, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(tracer_test.mSuite.runs[run_index].basePath,
                                    tracer_test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename + '.png', dpi = 300)
        plt.savefig(img_filename + '.pdf')
        plt.clf()
        tracer_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], tracer_test.outputPathBase):
    report_filename = os.path.join(tracer_test.outputPathBase,
                     "%s-report.%s" % (tracer_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(tracer_test, mResults, rGen, report_filename)
    html_filename = os.path.join(tracer_test.outputPathBase,
                     "%s-report.%s" % (tracer_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
