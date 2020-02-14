"""
1-D MINC column test
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
parser.add_argument("-np", type = int, default = 1, help = "number of processes")
parser.add_argument("-d", "--docker", action = "store_true",
                    help = "run via Docker (waiwera-dkr)")
args = parser.parse_args()
mpi = args.np > 1 and not args.docker
simulator = 'waiwera-dkr -np %d' % args.np if args.docker else 'waiwera'

model_name = 'minc_column'

AUTOUGH2_FIELDMAP = {
    'Vapour saturation': 'Vapour saturation',
    'Temperature': 'Temperature'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Enthalpy': 'source_enthalpy'}

model_dir = './run'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)

def minc_level_map(num_levels, num_minc_cells):
    # return mapping to reorder AUTOUGH2 cells into MINC levels, as in Waiwera
    m = list(range(geo.num_atmosphere_blocks, geo.num_blocks))
    for l in range(num_levels):
        level_map = list(range(geo.num_blocks + l,
                               geo.num_blocks + l + num_minc_cells,
                               num_levels))
        m += level_map
    return m

run_names = ['single', 'minc']

test_fields = ['Pressure', 'Temperature', 'Vapour saturation']
plot_fields = test_fields
test_source_fields = ["Enthalpy"]
field_scale = {'Pressure': 1.e5, 'Temperature': 1., 'Vapour saturation': 1., 'Enthalpy': 1.e3}
field_unit = {'Pressure': 'bar', 'Temperature': '$^{\circ}$C', 'Vapour saturation': '',
              'Enthalpy': 'kJ/kg'}

minc_column_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
minc_column_test.description = """1-D MINC column problem, production run starting
from steady-state solution, with single-porosity results included for comparison.
"""

for run_index, run_name in enumerate(run_names):
    base_path = os.path.realpath(model_dir)
    run_base_name = model_name + '_' + run_name
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                                fieldname_map = WAIWERA_FIELDMAP,
                                simulator = simulator,
                                basePath = base_path)
    model_run.jobParams['nproc'] = args.np
    minc_column_test.mSuite.addRun(model_run, run_name)

obs_cell_elev = -350.
obs_position = np.array([50., 50., obs_cell_elev])
obs_blk = geo.block_name_containing_point(obs_position)
obs_cell_index = geo.block_name_index[obs_blk] - geo.num_atmosphere_blocks
source_index = -1

minc_column_test.setupEmptyTestCompsList()
AUTOUGH2_result = {}

for run_index, run_name in enumerate(run_names):
    run_base_name = model_name + '_' + run_name
    results_filename = os.path.join(model_dir, run_base_name + ".listing")
    run_filename = run_base_name + '.json'
    inp = json.load(open(os.path.join(base_path, run_filename)))
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
    minc_column_test.addTestComp(run_index, "AUTOUGH2",
                      FieldWithinTolTC(fieldsToTest = test_fields,
                                       defFieldTol = 0.025,
                                       expected = AUTOUGH2_result[run_name],
                                       testOutputIndex = -1))

    minc_column_test.addTestComp(run_index, "AUTOUGH2 history",
                          HistoryWithinTolTC(fieldsToTest = test_fields,
                                             defFieldTol = 0.02,
                                             expected = AUTOUGH2_result[run_name],
                                             testCellIndex = obs_cell_index))

    minc_column_test.addTestComp(run_index, "AUTOUGH2 source",
                          HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                             defFieldTol = 0.01,
                                             expected = AUTOUGH2_result[run_name],
                                             testSourceIndex = source_index))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = minc_column_test.runTest(jrunner, createReports = True)

day = 24. * 60. * 60.

for field_name in plot_fields:
    scale = field_scale[field_name]
    unit = field_unit[field_name]
    for run_index, run_name in enumerate(run_names):
        title = run_name.replace('single', 'single porosity').replace('minc', 'MINC')
        t, var = minc_column_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t / day, var / scale, '-', label = 'Waiwera ' + title, zorder = 3)

        t, var = AUTOUGH2_result[run_name].getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t[::3] / day, var[::3] / scale, 's', label = 'AUTOUGH2 ' + title, zorder = 2)

        plt.xlabel('time (days)')
        plt.ylabel(field_name + ' (' + unit + ')')
        plt.legend(loc = 'best')
        plt.title(' '.join((field_name, 'history at production well')))
        img_filename_base = '_'.join((model_name, run_name, 'history', field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(minc_column_test.mSuite.runs[run_index].basePath,
                                    minc_column_test.mSuite.outputPathBase,
                                    img_filename_base)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    minc_column_test.mSuite.analysisImages.append(img_filename + '.png')

for field_name in test_source_fields:
    scale = field_scale[field_name]
    unit = field_unit[field_name]
    for run_index, run_name in enumerate(run_names):
        title = run_name.replace('single', 'single porosity').replace('minc', 'MINC')
        t, var = minc_column_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t / day, var / scale, '-', label = 'Waiwera ' + title, zorder = 3)

        t, var = AUTOUGH2_result[run_name].getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t[::3] / day, var[::3] / scale, 's', label = 'AUTOUGH2 ' + title, zorder = 2)

        plt.xlabel('time (days)')
        plt.ylabel(field_name + ' (' + unit + ')')
        plt.legend(loc = 'best')
        plt.title(' '.join((field_name, 'history at production well')))
        img_filename_base = '_'.join((model_name, run_name, 'history', field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(minc_column_test.mSuite.runs[run_index].basePath,
                                    minc_column_test.mSuite.outputPathBase,
                                    img_filename_base)
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    minc_column_test.mSuite.analysisImages.append(img_filename + '.png')

z = [lay.centre for lay in geo.layerlist[geo.num_atmosphere_blocks:]]
n = len(z)

for field_name in plot_fields:
    scale = field_scale[field_name]
    unit = field_unit[field_name]
    for run_index, run_name in enumerate(run_names):
        result = minc_column_test.mSuite.resultsList[run_index]
        var = result.getFieldAtOutputIndex(field_name, -1) / scale
        title = run_name.replace('single', 'single porosity').replace('minc', 'MINC')
        plt.plot(var[:n], z, '-', label = 'Waiwera ' + title)
        var = AUTOUGH2_result[run_name].getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(var[:n], z, 's', label = 'AUTOUGH2 ' + title)
    plt.xlabel(field_name + ' (' + unit + ')')
    plt.ylabel('z (m)')
    plt.title(' '. join(['Final', field_name.lower(), 'profile']))
    img_filename_base = '_'.join((model_name, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(minc_column_test.mSuite.runs[0].basePath,
                                minc_column_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.legend(loc = 'best')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    minc_column_test.mSuite.analysisImages.append(img_filename + '.png')

# generate report:

for rGen in getGenerators(["RST"], minc_column_test.outputPathBase):
    report_filename = os.path.join(minc_column_test.outputPathBase,
                     "%s-report.%s" % (minc_column_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(minc_column_test, mResults, rGen, report_filename)
    html_filename = os.path.join(minc_column_test.outputPathBase,
                     "%s-report.%s" % (minc_column_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
