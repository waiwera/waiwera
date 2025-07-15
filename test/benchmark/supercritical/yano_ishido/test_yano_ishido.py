"""
Yano-Ishido supercritical radial problem - temperature cases
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest

from credo.jobrunner import SimpleJobRunner
from credo.waiwera import WaiweraModelRun
from credo.modelresult import HistoryDataResult

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from credo.systest import HistoryWithinTolTC

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'

import numpy as np
from docutils.core import publish_file

parser = argparse.ArgumentParser()
parser.add_argument("-np", type = int, default = 1, help = "number of processes")
parser.add_argument("-d", "--docker", action = "store_true",
                    help = "run via Docker (waiwera-dkr)")
args = parser.parse_args()
mpi = args.np > 1 and not args.docker
simulator = 'waiwera-dkr -np %d' % args.np if args.docker else 'waiwera'

model_name = 'yanoish'

WAIWERA_FIELDMAP = {'Pressure': 'fluid_pressure'}

model_dir = './run'
data_dir = './data'

temps = np.arange(200., 550., 50.)
run_names = [str(int(T0)) for T0 in temps]
run_symbols = ['.', '^', 's', 'v', 'D', '+', 'o']
run_colours = ['r', 'g', 'm', 'k', 'y', 'c', 'b']
obs_cell_index = 0
hrs = 60. * 60.

test_fields = ['Pressure']
field_scale = {'Pressure': 1.e6}
field_unit = {'Pressure': 'MPa'}

yano_ishido_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
yano_ishido_test.description = """ Yano-Ishido supercritical radial problem, from Yano and Ishido
(1998), results digitised from figure 8 """

for run_index, run_name in enumerate(run_names):
    run_base_name = model_name + '_' + run_name
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                                fieldname_map = WAIWERA_FIELDMAP,
                                simulator = simulator,
                                basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = args.np
    yano_ishido_test.mSuite.addRun(model_run, run_name)

yano_ishido_test.setupEmptyTestCompsList()

digitised_result = {}

for run_index, run_name in enumerate(run_names):
    data = {}
    for field_name in test_fields:
        data_filename = '_'.join(('fig8', run_name))
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        data[field_name, obs_cell_index] = np.loadtxt(data_filename)
        data[field_name, obs_cell_index][:,0] *= hrs
        data[field_name, obs_cell_index][:,1] *= field_scale[field_name]
        digitised_result[run_name] = HistoryDataResult(run_name, data)
        yano_ishido_test.addTestComp(run_index, 'T0 = ' + run_name + ' deg C',
                              HistoryWithinTolTC(fieldsToTest = test_fields,
                                                 defFieldTol = 0.03,
                                                 expected = digitised_result[run_name],
                                                 testCellIndex = obs_cell_index,
                                                 orthogonalError = True,
                                                 logx = False))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = yano_ishido_test.runTest(jrunner, createReports = True)

symbols = dict(zip(run_names, run_symbols))
colours = dict(zip(run_names, run_colours))
tmin = 500
for run_index, run_name in enumerate(run_names):
    for field_name in test_fields:

        t, var = yano_ishido_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtCell(field_name, obs_cell_index)
        ix = np.where(t > tmin)
        t, var = t[ix], var[ix]
        plt.semilogx(t / hrs, var / field_scale[field_name],
                     symbols[run_name], color = colours[run_name], zorder = 4)

        result = digitised_result[run_name]
        t, var = result.getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.semilogx(t / hrs, var / field_scale[field_name],
                     '-', color = colours[run_name],
                     label = run_name + '$^{\circ}$C')

plt.xlabel('time (hr)')
plt.ylabel(field_name + ' (' + field_unit[field_name] + ')')

plt.legend(loc = 'best')
plt.title('Yano-Ishido pressure history')
img_filename_base = '_'.join((model_name, field_name))
img_filename_base = img_filename_base.replace(' ', '_')
img_filename = os.path.join(yano_ishido_test.mSuite.runs[run_index].basePath,
                            yano_ishido_test.mSuite.outputPathBase,
                            img_filename_base)
plt.tight_layout(pad = 3.)
plt.savefig(img_filename + '.png', dpi = 300)
plt.savefig(img_filename + '.pdf')
plt.clf()
yano_ishido_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], yano_ishido_test.outputPathBase):
    report_filename = os.path.join(yano_ishido_test.outputPathBase,
                     "%s-report.%s" % (yano_ishido_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(yano_ishido_test, mResults, rGen, report_filename)
    html_filename = os.path.join(yano_ishido_test.outputPathBase,
                     "%s-report.%s" % (yano_ishido_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
