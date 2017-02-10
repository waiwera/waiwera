"""
Model Intercomparison Study problem 2
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
from t2thermo import *
from scipy.special import expi

from docutils.core import publish_file

P0 = 90.e5
T0 = 260.
qm = -14.
k = 1.e-14
phi = 0.2
h = 100.
ps = sat(T0)
pav = 81.e5
rhol, ul = cowat(T0, pav)
vl = visw(T0, pav, ps) / rhol
comp = 1.5805e-9
D = k / (vl * phi * rhol * comp)
def theis(pos, t):
    r = pos[0]
    tol = 1.e-6
    if t < tol: return P0
    else:
        s = -r * r / (4. * D * t)
        d = 4. * np.pi * k * h / vl
        return P0 - qm * expi(s) / d

model_name = 'problem2'

WAIWERA_FIELDMAP = {'Pressure': 'fluid_pressure'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)

num_procs = 1
run_name = 'case a'
run_index = 0
test_fields = ["Pressure"]

problem2_test = SciBenchmarkTest(model_name + "_test", nproc = num_procs)
problem2_test.description = """Model Intercomparison Study problem 2
(Theis problem)"""

run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                          fieldname_map = WAIWERA_FIELDMAP,
                          simulator = 'waiwera',
                          basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = num_procs
problem2_test.mSuite.addRun(model_run, run_name)
problem2_test.setupEmptyTestCompsList()

run_base_name = model_name
results_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", results_filename,
                                 geo_filename = t2geo_filename)

obs_radii = [0.15, 0.5, 1.]
for r in obs_radii:
    obspt = 'r = ' + str(r) + ' m'
    pos = np.array([r, 0., -0.5 * h])
    obs_blk = geo.block_name_containing_point(pos)
    obs_cell_index = geo.block_name_index[obs_blk]
    problem2_test.addTestComp(run_index, "AUTOUGH2 " + obspt,
                              HistoryWithinTolTC(fieldsToTest = test_fields,
                                                 defFieldTol = 1.e-5,
                                                 expected = AUTOUGH2_result,
                                                 testCellIndex = obs_cell_index))
    problem2_test.addTestComp(run_index, "Analytical solution " + obspt,
                              HistoryWithinTolTC(fieldsToTest = test_fields,
                                                 defFieldTol = 1.e-2,
                                                 expected = theis,
                                                 testCellIndex = obs_cell_index))

jrunner = SimpleJobRunner(mpi = True)
testResult, mResults = problem2_test.runTest(jrunner, createReports = True)

# plots:
field_name = "Pressure"
bar, unit = 1.e5, 'bar'
for r in obs_radii:
    obspt = 'r = ' + str(r) + ' m'
    pos = np.array([r, 0., -0.5 * h])
    obs_blk = geo.block_name_containing_point(pos)
    obs_cell_index = geo.block_name_index[obs_blk]
    t = problem2_test.testComps[run_index]["AUTOUGH2 " + obspt].times
    var = problem2_test.mSuite.resultsList[run_index].getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.semilogx(t, var / bar, '-', label = 'Waiwera')
    t = AUTOUGH2_result.getTimes()
    var = AUTOUGH2_result.getFieldHistoryAtCell(field_name, obs_cell_index)
    plt.semilogx(t, var / bar, '+', label = 'AUTOUGH2')
    var = np.array([theis(pos, ti) for ti in t])
    plt.semilogx(t, var / bar, ':', label = 'analytical')
    plt.xlabel('time (s)')
    plt.ylabel(field_name + '(' + unit + ')')
    plt.legend()
    plt.title(' '.join((model_name, run_name, field_name.lower(),
                        'results:', obspt)))
    img_filename_base = '_'.join((model_name, run_name, obspt, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(problem2_test.mSuite.runs[run_index].basePath,
                                problem2_test.mSuite.outputPathBase,
                                img_filename_base + '.png')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename)
    plt.clf()
    problem2_test.mSuite.analysisImages.append(img_filename)
    
# generate report:

for rGen in getGenerators(["RST"], problem2_test.outputPathBase):
    report_filename = os.path.join(problem2_test.outputPathBase,
                     "%s-report.%s" % (problem2_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(problem2_test, mResults, rGen, report_filename)
    html_filename = os.path.join(problem2_test.outputPathBase,
                     "%s-report.%s" % (problem2_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")

