"""
Model Intercomparison Study problem 2
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest
from similarity import DigitisedSimilarityResult, \
    SimilaritySolutionWithinTolTC

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import ModelResult
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from mulgrids import mulgrid

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'

import numpy as np
from t2thermo import *
from scipy.special import expi

from docutils.core import publish_file

parser = argparse.ArgumentParser()
parser.add_argument("-np", type = int, default = 1, help = "number of processes")
parser.add_argument("-d", "--docker", action = "store_true",
                    help = "run via Docker (waiwera-dkr)")
args = parser.parse_args()
mpi = args.np > 1 and not args.docker
simulator = 'waiwera-dkr -np %d' % args.np if args.docker else 'waiwera'

model_name = 'problem2'

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Liquid saturation': 'fluid_liquid_saturation'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)

run_names = ['a', 'b', 'c']
test_fields = {
    'a': ["Pressure"],
    'b': ["Pressure", "Liquid saturation"],
    'c': ["Pressure", "Liquid saturation"]}
field_scale = {"Pressure": 1.e-5, "Liquid saturation": 1.}
field_unit = {"Pressure": "bar", "Liquid saturation": ""}
obs_radii = [0.5, 1.]
obs_cell_indices = []
thickness = 100.
for r in obs_radii:
    pos = np.array([r, 0., -0.5 * thickness])
    obs_blk = geo.block_name_containing_point(pos)
    obs_cell_indices.append(geo.block_name_index[obs_blk])

# Theis solution for case a:
day = 24. * 60. * 60.
P0 = 90.e5
T0 = 260.
qm = -14.
k = 1.e-14
phi = 0.2
ps = sat(T0)
pav = 81.e5
rhol, ul = cowat(T0, pav)
vl = visw(T0, pav, ps) / rhol
comp = 1.5805e-9
D = k / (vl * phi * rhol * comp)
def theis_sim(sim):
    tol = 1.e-5
    if sim < tol: return P0
    else:
        s = 1. / (4. * D * sim)
        d = 4. * np.pi * k * thickness / vl
        return P0 - qm * expi(s) / d
def theis(pos, t):
    r = pos[0]
    return theis_sim(t / (r * r))

expected = {
    'a': theis,
    'b': DigitisedSimilarityResult('semi-analytical',
                                   {'Pressure': 'data/problem2b_pressure.dat',
                                    'Liquid saturation': 'data/problem2b_saturation.dat'}),
    'c': DigitisedSimilarityResult('reference',
                                   {'Pressure': 'data/problem2c_pressure_s-cubed.dat',
                                    'Liquid saturation': 'data/problem2c_saturation_s-cubed.dat'})
}
ref_result_name = {'a': 'Theis solution',
                   'b': 'semi-analytical',
                   'c': 'S-Cubed'}
symbol = {'Theis solution': 'k--',
          'semi-analytical': 'k--',
          'S-Cubed': 'o'}
AUTOUGH2_field_tol = {'a': 1.e-4, 'b': 1.e-4, 'c': 1.e-2}
ref_field_tol = {'a': 2.e-2, 'b': 2.e-2, 'c': 3.5e-2}
case_c_semi_analytical = DigitisedSimilarityResult('semi-analytical',
                                   {'Pressure': 'data/problem2c_pressure_analytical.dat',
                                    'Liquid saturation': 'data/problem2c_saturation_analytical.dat'})
case_b_s_cubed = DigitisedSimilarityResult('S-Cubed',
                                   {'Pressure': 'data/problem2b_pressure_s-cubed.dat',
                                    'Liquid saturation': 'data/problem2b_saturation_s-cubed.dat'})

# NB here I have used the semi-analytical results as reference for
# case b, as they are tolerably close to all the simulation results,
# but for case c they are not, so I have used the S-Cubed results as
# reference there. However both semi-analytical and S-Cubed results
# are shown on plots for cases b and c.

problem2_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
problem2_test.description = """Model Intercomparison Study problem 2
(case a: Theis problem, case b: radial two-phase production, case c: radial flashing front)"""

for run_index, run_name in enumerate(run_names):

    run_base_name = model_name + run_name
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                              fieldname_map = WAIWERA_FIELDMAP,
                              simulator = simulator,
                              basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = args.np
    problem2_test.mSuite.addRun(model_run, run_name)

problem2_test.setupEmptyTestCompsList()
AUTOUGH2_result = {}

for run_index, run_name in enumerate(run_names):

    run_base_name = model_name + run_name
    results_filename = os.path.join(model_dir, run_base_name + ".listing")
    AUTOUGH2_result[run_name] = T2ModelResult("AUTOUGH2", results_filename,
                                    geo_filename = t2geo_filename)

    problem2_test.addTestComp(run_index, "AUTOUGH2",
                              SimilaritySolutionWithinTolTC(
                                  fieldsToTest = test_fields[run_name],
                                  defFieldTol = AUTOUGH2_field_tol[run_name],
                                  expected = AUTOUGH2_result[run_name],
                                  testCellIndices = obs_cell_indices))
    problem2_test.addTestComp(run_index,
                              ref_result_name[run_name],
                              SimilaritySolutionWithinTolTC(
                                  fieldsToTest = test_fields[run_name],
                                  defFieldTol = ref_field_tol[run_name],
                                  fieldTols = {"Liquid saturation": 6.e-2},
                                  absoluteErrorTol = 0.01,
                                  expected = expected[run_name],
                                  testCellIndices = obs_cell_indices))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = problem2_test.runTest(jrunner, createReports = True)

# plots:
for run_index, run_name in enumerate(run_names):
    for field_name in test_fields[run_name]:
        scale = field_scale[field_name]
        sims = []
        var_waiwera, var_AUTOUGH2 = [], []
        for r in obs_radii:
            r2 = r * r
            obspt = 'r = ' + str(r) + ' m'
            pos = np.array([r, 0., -0.5 * thickness])
            obs_blk = geo.block_name_containing_point(pos)
            obs_cell_index = geo.block_name_index[obs_blk]
            t, var = problem2_test.mSuite.resultsList[run_index].\
                  getFieldHistoryAtCell(field_name, obs_cell_index)
            sim = t / r2
            sims += list(sim)
            var_waiwera += list(var)
            t, var = AUTOUGH2_result[run_name].getFieldHistoryAtCell(field_name,
                                                                  obs_cell_index)
            var_AUTOUGH2 += list(var)
        iord = np.argsort(np.array(sims))
        plt.semilogx(np.array(sims)[iord] / day, np.array(var_waiwera)[iord] * scale,
                     '-', label = 'Waiwera', zorder = 4)
        plt.semilogx(np.array(sims) / day, np.array(var_AUTOUGH2) * scale,
                     's', label = 'AUTOUGH2', zorder = 3)
        if run_name == 'a':
            sims = np.logspace(-5, 1.8, 20) * day
            var = np.array([theis_sim(sim) for sim in sims])
        else:
            sims = expected[run_name].getSimilarityVariables(field_name)
            var = expected[run_name].getSimilarityValues(field_name)
        ref_name = ref_result_name[run_name]
        plt.semilogx(sims / day, var * scale, symbol[ref_name], label = ref_name)
        if run_name == 'b':
            sims = case_b_s_cubed.getSimilarityVariables(field_name)
            var = case_b_s_cubed.getSimilarityValues(field_name)
            plt.semilogx(sims / day, var * scale, 'o', label = 'S-Cubed')
        elif run_name == 'c':
            sims = case_c_semi_analytical.getSimilarityVariables(field_name)
            var = case_c_semi_analytical.getSimilarityValues(field_name)
            plt.semilogx(sims / day, var * scale, 'k--', label = 'semi-analytical')
        plt.xlabel('t /$r^2$ (day/$m^2$)')
        plt.ylabel(field_name + ' (' + field_unit[field_name] + ')')
        plt.legend(loc = 'best')
        plt.title(' '.join((model_name + run_name, field_name.lower(),
                            'results')))
        img_filename_base = '_'.join((model_name, run_name, obspt, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(problem2_test.mSuite.runs[run_index].basePath,
                                    problem2_test.mSuite.outputPathBase,
                                    img_filename_base)
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename +'.png', dpi = 300)
        plt.savefig(img_filename +'.pdf')
        plt.clf()
        problem2_test.mSuite.analysisImages.append(img_filename + '.png')
    
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

