"""
3-D MINC production problem with well on deliverability, and
time-dependent productivity index
"""

import argparse
import os
import sys

from credo.systest import SciBenchmarkTest

from credo.jobrunner import SimpleJobRunner
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from credo.systest import FieldWithinTolTC
from credo.systest import HistoryWithinTolTC

from mulgrids import mulgrid
from t2listing import t2listing

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'

import numpy as np
import json
from docutils.core import publish_file

parser = argparse.ArgumentParser()
parser.add_argument("-np", help = "number of processes")
args = parser.parse_args()
if args.np: num_procs = int(args.np)
else: num_procs = 1

model_name = 'minc_3d'

AUTOUGH2_FIELDMAP = {
    'Vapour saturation': 'Vapour saturation',
    'Temperature': 'Temperature'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Enthalpy': 'source_enthalpy',
    'Generation rate': 'source_rate'}

model_dir = './run'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)

def minc_level_map(num_levels, num_minc_cells):
    # return mapping to reorder AUTOUGH2 cells into MINC levels, as in Waiwera
    m = range(geo.num_atmosphere_blocks, geo.num_blocks)
    for l in range(num_levels):
        level_map = list(range(geo.num_blocks + l,
                               geo.num_blocks + l + num_minc_cells,
                               num_levels))
        m += level_map
    return m

run_names = ['MINC']

test_fields = ['Pressure', 'Temperature', 'Vapour saturation']
plot_fields = test_fields
test_source_fields = ['Enthalpy', 'Generation rate']
field_scale = {'Pressure': 1.e5, 'Temperature': 1., 'Vapour saturation': 1.,
               'Enthalpy': 1.e3, 'Generation rate': 1.}
field_unit = {'Pressure': 'bar', 'Temperature': '$^{\circ}$C', 'Vapour saturation': '',
              'Enthalpy': 'kJ/kg', 'Generation rate': 'kg/s'}

minc_production_test = SciBenchmarkTest(model_name + "_test", nproc = num_procs)
minc_production_test.description = """3-D MINC production problem, with well on deliverability,
time-dependent productivity index and steam flow limiter.
"""

for run_index, run_name in enumerate(run_names):
    base_path = os.path.realpath(model_dir)
    run_base_name = model_name
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                                fieldname_map = WAIWERA_FIELDMAP,
                                simulator = 'waiwera',
                                basePath = base_path)
    model_run.jobParams['nproc'] = num_procs
    minc_production_test.mSuite.addRun(model_run, run_name)

obs_cell_elev = -1000.
obs_position = np.array([10., 10., obs_cell_elev])
obs_blk = geo.block_name_containing_point(obs_position)
obs_cell_index = geo.block_name_index[obs_blk] - geo.num_atmosphere_blocks
source_index = -1

minc_production_test.setupEmptyTestCompsList()
AUTOUGH2_result = {}

for run_index, run_name in enumerate(run_names):
    run_base_name = model_name
    results_filename = os.path.join(model_dir, run_base_name + ".listing")
    run_filename = run_base_name + '.json'
    inp = json.load(file(os.path.join(base_path, run_filename)))
    if 'minc' in inp['mesh']:
        num_levels = len(inp['mesh']['minc']['geometry']['matrix']['volume'])
    else: num_levels = 0
    lst = t2listing(results_filename)
    lst.last()
    num_minc_cells = lst.element.num_rows - geo.num_blocks
    AUTOUGH2_result[run_name] = T2ModelResult("AUTOUGH2", results_filename,
                                              geo_filename = t2geo_filename,
                                              fieldname_map = AUTOUGH2_FIELDMAP,
                                              ordering_map = minc_level_map(num_levels,
                                                                            num_minc_cells))
    minc_production_test.addTestComp(run_index, "AUTOUGH2",
                      FieldWithinTolTC(fieldsToTest = test_fields,
                                       defFieldTol = 0.01,
                                       expected = AUTOUGH2_result[run_name],
                                       testOutputIndex = -1))

    minc_production_test.addTestComp(run_index, "AUTOUGH2 history",
                          HistoryWithinTolTC(fieldsToTest = test_fields,
                                             defFieldTol = 0.01,
                                             expected = AUTOUGH2_result[run_name],
                                             testCellIndex = obs_cell_index))

    minc_production_test.addTestComp(run_index, "AUTOUGH2 source",
                          HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                             defFieldTol = 0.01,
                                             expected = AUTOUGH2_result[run_name],
                                             testSourceIndex = source_index))

jrunner = SimpleJobRunner(mpi = True)
testResult, mResults = minc_production_test.runTest(jrunner, createReports = True)

yr = 365. * 24. * 60. * 60.

for field_name in plot_fields:
    scale = field_scale[field_name]
    unit = field_unit[field_name]
    for run_index, run_name in enumerate(run_names):
        t, var = minc_production_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t / yr, var / scale, '-', label = 'Waiwera', zorder = 3)

        t, var = AUTOUGH2_result[run_name].getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t[::3] / yr, var[::3] / scale, 's', label = 'AUTOUGH2', zorder = 2)

        plt.xlabel('time (years)')
        plt.ylabel(field_name + ' (' + unit + ')')
        plt.legend(loc = 'best')
        plt.title(' '.join((field_name, 'history at production well')))
        img_filename_base = '_'.join((model_name, run_name, 'history', field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(minc_production_test.mSuite.runs[run_index].basePath,
                                    minc_production_test.mSuite.outputPathBase,
                                    img_filename_base)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    minc_production_test.mSuite.analysisImages.append(img_filename + '.png')

for field_name in test_source_fields:
    scale = field_scale[field_name]
    unit = field_unit[field_name]
    for run_index, run_name in enumerate(run_names):
        t, var = minc_production_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t / yr, var / scale, '-', label = 'Waiwera', zorder = 3)

        t, var = AUTOUGH2_result[run_name].getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t[::3] / yr, var[::3] / scale, 's', label = 'AUTOUGH2', zorder = 2)

        plt.xlabel('time (years)')
        plt.ylabel(field_name + ' (' + unit + ')')
        plt.legend(loc = 'best')
        plt.title(' '.join((field_name, 'history at production well')))
        img_filename_base = '_'.join((model_name, run_name, 'history', field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(minc_production_test.mSuite.runs[run_index].basePath,
                                    minc_production_test.mSuite.outputPathBase,
                                    img_filename_base)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    minc_production_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], minc_production_test.outputPathBase):
    report_filename = os.path.join(minc_production_test.outputPathBase,
                     "%s-report.%s" % (minc_production_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(minc_production_test, mResults, rGen, report_filename)
    html_filename = os.path.join(minc_production_test.outputPathBase,
                     "%s-report.%s" % (minc_production_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
