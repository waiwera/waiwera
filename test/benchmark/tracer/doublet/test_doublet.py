"""
1-D single-phase doublet problem with tracer injection and production
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

model_name = 'doublet'
    
AUTOUGH2_FIELDMAP = {
    'Pressure': 'Pressure',
    'Tracer mass fraction': 'Tracer/liquid'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Tracer mass fraction': 'tracer_1'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
map_out_bdy = list(range(0, geo.num_blocks))

run_index = 0
run_name = 'run'
run_base_name = model_name
output_indices = [0, 1, 2, 3, 6, 8]

test_fields = ['Tracer mass fraction']
field_unit = {'Tracer mass fraction': '10$^{-6}$'}
field_scale = {'Tracer mass fraction': 1.e-6}

doublet_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
doublet_test.description = """1-D single-phase doublet problem with tracer injection
and production
"""

run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                            fieldname_map = WAIWERA_FIELDMAP,
                            simulator = simulator,
                            basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = args.np
doublet_test.mSuite.addRun(model_run, run_name)

doublet_test.setupEmptyTestCompsList()

results_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", results_filename,
                                geo_filename = t2geo_filename,
                                fieldname_map = AUTOUGH2_FIELDMAP,
                                ordering_map = map_out_bdy)
for output_index in output_indices:
    doublet_test.addTestComp(run_index, "AUTOUGH2_%d" % output_index,
                             FieldWithinTolTC(fieldsToTest = test_fields,
                                              defFieldTol = 1.e-3,
                                              expected = AUTOUGH2_result,
                                              testOutputIndex = output_index))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = doublet_test.runTest(jrunner, createReports = True)

result = doublet_test.mSuite.resultsList[run_index]
t = result.getTimes()
x = np.array([col.centre[0] for col in geo.columnlist])

for field_name in test_fields:
    for output_index in output_indices:
        output_str = str(output_index)
        var = result.getFieldAtOutputIndex(field_name,
                                           output_index) / field_scale[field_name]
        plt.plot(x, var, 'b-', label = 'Waiwera', zorder = 2)
        var = AUTOUGH2_result.getFieldAtOutputIndex(field_name,
                                                    output_index) / field_scale[field_name]
        plt.plot(x, var, 'gs', label = 'AUTOUGH2', zorder = 1)
        plt.xlabel('x (m)')
        plt.ylabel('%s (%s)' % (field_name, field_unit[field_name]))
        img_filename_base = '_'.join((model_name, run_name, field_name, output_str))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(doublet_test.mSuite.runs[run_index].basePath,
                                    doublet_test.mSuite.outputPathBase,
                                    img_filename_base)
        plt.legend(loc = 'best')
        plt.title('%s at time %2.0f days' % (field_name, t[output_index] / day))
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename + '.png', dpi = 300)
        plt.savefig(img_filename + '.pdf')
        plt.clf()
        doublet_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], doublet_test.outputPathBase):
    report_filename = os.path.join(doublet_test.outputPathBase,
                     "%s-report.%s" % (doublet_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(doublet_test, mResults, rGen, report_filename)
    html_filename = os.path.join(doublet_test.outputPathBase,
                     "%s-report.%s" % (doublet_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
