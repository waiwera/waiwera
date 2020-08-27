"""
One-cell tracer decay
"""

import argparse
import os
import sys
from math import exp
from functools import partial

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest

from credo.jobrunner import SimpleJobRunner
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from credo.systest import HistoryWithinTolTC

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

model_name = 'decay'
tracers = ['no_decay', 'constant', 'temperature']

T = 60.
X = 0.001
k0 = 1.e-6
Ea = 2.e3
R = 8.3144598
Tk = T + 273.15
decay_rates = [0, k0, k0 * exp(-Ea / (R * Tk))]

def exact(pos, t, k): return X * exp(-k * t)

model_dir = './run'
data_dir = './data'

run_name = 'run'
run_index = 0

decay_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
decay_test.description = """One-cell tracer decay
"""
run_filename = model_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                            simulator = simulator,
                            basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = args.np
decay_test.mSuite.addRun(model_run, run_name)

decay_test.setupEmptyTestCompsList()

for tracer, decay_rate in zip(tracers, decay_rates):
    decay_test.addTestComp(run_index, tracer,
                           HistoryWithinTolTC(fieldsToTest = ['tracer_' + tracer],
                                              defFieldTol = 1.e-2,
                                              absoluteErrorTol = 1.e-4,
                                              expected = partial(exact, k = decay_rate),
                                              testCellIndex = 0))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = decay_test.runTest(jrunner, createReports = True)

for tracer, decay_rate in zip(tracers, decay_rates):
    result = decay_test.mSuite.resultsList[run_index]
    t, x = decay_test.mSuite.resultsList[run_index].\
             getFieldHistoryAtCell('tracer_' + tracer, 0)
    plt.plot(t / day, x, 'b-', label = 'Waiwera', zorder = 2)
    xexact = [exact(0, ti, decay_rate) for ti in t]
    plt.plot(t / day, xexact, 'gs', label = 'exact', zorder = 1)
    plt.xlabel('t (day)')
    plt.ylabel('tracer mass fraction')
    img_filename_base = '_'.join((model_name, tracer))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(decay_test.mSuite.runs[run_index].basePath,
                                decay_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.legend(loc = 'best')
    plt.title(tracer)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    decay_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], decay_test.outputPathBase):
    report_filename = os.path.join(decay_test.outputPathBase,
                     "%s-report.%s" % (decay_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(decay_test, mResults, rGen, report_filename)
    html_filename = os.path.join(decay_test.outputPathBase,
                     "%s-report.%s" % (decay_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
