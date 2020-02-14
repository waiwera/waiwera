"""
1-D MINC doublet test
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

from mulgrids import mulgrid

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

model_name = 'minc_1d'

AUTOUGH2_FIELDMAP = {
    'Vapour saturation': 'Vapour saturation',
    'Temperature': 'Temperature'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
num_cols = geo.num_columns

def minc_level_map(num_levels):
    # return mapping to reorder AUTOUGH2 cells into MINC levels, as in Waiwera
    m = list(range(num_cols))
    for l in range(num_levels):
        level_map = list(range(num_cols + l,
                               num_cols + l + num_levels * num_cols,
                               num_levels))
        m += level_map
    return m

spacings = [50, 100, 200] # fracture spacings
run_names = ['single'] + [str(s) for s in spacings]

test_fields = ['Pressure', 'Temperature', 'Vapour saturation']
plot_fields = test_fields
field_scale = {'Pressure': 1.e5, 'Temperature': 1., 'Vapour saturation': 1.}
field_unit = {'Pressure': 'bar', 'Temperature': '$^{\circ}$C', 'Vapour saturation': ''}

minc_1d_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
minc_1d_test.description = """1-D MINC doublet problem
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
    minc_1d_test.mSuite.addRun(model_run, run_name)

minc_1d_test.setupEmptyTestCompsList()
AUTOUGH2_result = {}

for run_index, run_name in enumerate(run_names):
    run_base_name = model_name + '_' + run_name
    results_filename = os.path.join(model_dir, run_base_name + ".listing")
    run_filename = run_base_name + '.json'
    inp = json.load(open(os.path.join(base_path, run_filename)))
    if 'minc' in inp['mesh']:
        num_levels = len(inp['mesh']['minc']['geometry']['matrix']['volume'])
    else: num_levels = 0
    AUTOUGH2_result[run_name] = T2ModelResult("AUTOUGH2", results_filename,
                                              geo_filename = t2geo_filename,
                                              fieldname_map = AUTOUGH2_FIELDMAP,
                                              ordering_map = minc_level_map(num_levels))
    minc_1d_test.addTestComp(run_index, "AUTOUGH2",
                      FieldWithinTolTC(fieldsToTest = test_fields,
                                       defFieldTol = 2.e-3,
                                       expected = AUTOUGH2_result[run_name],
                                       testOutputIndex = -1))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = minc_1d_test.runTest(jrunner, createReports = True)

x = [col.centre[0] for col in geo.columnlist]
n = geo.num_columns

for field_name in plot_fields:
    scale = field_scale[field_name]
    unit = field_unit[field_name]
    for run_index, run_name in enumerate(run_names):
        result = minc_1d_test.mSuite.resultsList[run_index]
        var = result.getFieldAtOutputIndex(field_name, -1) / scale
        if run_index == 0: title = 'Single porosity'
        else: title = 'MINC (' + run_name + ' m)'
        plt.plot(x, var[:n], '-', label = 'Waiwera ' + title)
        var = AUTOUGH2_result[run_name].getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(x, var[:n], 's', label = 'AUTOUGH2 ' + title)
    plt.xlabel('x (m)')
    plt.ylabel(field_name + ' (' + unit + ')')
    plt.title(field_name)
    img_filename_base = '_'.join((model_name, field_name))
    img_filename_base = img_filename_base.replace(' ', '_')
    img_filename = os.path.join(minc_1d_test.mSuite.runs[0].basePath,
                                minc_1d_test.mSuite.outputPathBase,
                                img_filename_base)
    plt.legend(loc = 'best')
    plt.tight_layout(pad = 3.)
    plt.savefig(img_filename + '.png', dpi = 300)
    plt.savefig(img_filename + '.pdf')
    plt.clf()
    minc_1d_test.mSuite.analysisImages.append(img_filename + '.png')
        
# generate report:

for rGen in getGenerators(["RST"], minc_1d_test.outputPathBase):
    report_filename = os.path.join(minc_1d_test.outputPathBase,
                     "%s-report.%s" % (minc_1d_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(minc_1d_test, mResults, rGen, report_filename)
    html_filename = os.path.join(minc_1d_test.outputPathBase,
                     "%s-report.%s" % (minc_1d_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
