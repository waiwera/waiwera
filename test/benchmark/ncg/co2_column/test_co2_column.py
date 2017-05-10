"""
CO2 column test from O'Sullivan et al. (1985)
"""

import os

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
from docutils.core import publish_file

def total_CO2_mass_fraction(mResult, index):
    """Returns total CO2 mass fraction across all phases, from Waiwera
    results."""
    sl = mResult.getFieldAtOutputIndex('fluid_liquid_saturation', index)
    sv = mResult.getFieldAtOutputIndex('fluid_vapour_saturation', index)
    rhol = mResult.getFieldAtOutputIndex('fluid_liquid_density', index)
    rhov = mResult.getFieldAtOutputIndex('fluid_vapour_density', index)
    xl = mResult.getFieldAtOutputIndex('fluid_liquid_CO2_mass_fraction', index)
    xv = mResult.getFieldAtOutputIndex('fluid_vapour_CO2_mass_fraction', index)
    return (sl * rhol * xl + sv * rhov * xv) / (sl * rhol + sv * rhov)
    
model_name = 'co2_column'

AUTOUGH2_FIELDMAP = {
    'Vapour saturation': 'Gas saturatio',
    'CO2 mass fraction': 'CO2 mass fractio'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'CO2 mass fraction': total_CO2_mass_fraction}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
map_out_atm = range(geo.num_atmosphere_blocks, geo.num_blocks)

num_procs = 1
CO2_mass_fractions = [0, 0.1, 1, 5] # percent
run_names = [str(xgp) for xgp in CO2_mass_fractions]
obs_cell_index = 0

test_fields = ['Pressure', 'Temperature', 'Vapour saturation', 'CO2 mass fraction']
plot_fields = ['Vapour saturation', 'CO2 mass fraction']
field_scale = {'Temperature': 1., 'Vapour saturation': 1., 'CO2 mass fraction': 1.}
field_unit = {'Temperature': '$^{\circ}$C', 'Vapour saturation': '', 'CO2 mass fraction': ''}

co2_column_test = SciBenchmarkTest(model_name + "_test", nproc = num_procs)
co2_column_test.description = """Vertical column CO2 test from O'Sullivan et al. (1985). The original
problem (figs 10 and 11) specified CO2 input in terms of partial pressure, though the text indicated
that mass fraction (%) may actually have been used. Here, CO2 input is in terms of mass fraction. The
results, even for AUTOUGH2, do not match those in the original paper; however, the problem is still
useful as a benchmark of Waiwera against AUTOUGH2.
"""

for run_index, run_name in enumerate(run_names):
    run_base_name = model_name + '_' + run_name
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                                fieldname_map = WAIWERA_FIELDMAP,
                                simulator = 'waiwera',
                                basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = num_procs
    co2_column_test.mSuite.addRun(model_run, run_name)

co2_column_test.setupEmptyTestCompsList()
AUTOUGH2_result = {}

for run_index, run_name in enumerate(run_names):
    run_base_name = model_name + '_' + run_name
    results_filename = os.path.join(model_dir, run_base_name + ".listing")
    AUTOUGH2_result[run_name] = T2ModelResult("AUTOUGH2", results_filename,
                                              geo_filename = t2geo_filename,
                                              fieldname_map = AUTOUGH2_FIELDMAP,
                                              ordering_map = map_out_atm)
    co2_column_test.addTestComp(run_index, "AUTOUGH2",
                      FieldWithinTolTC(fieldsToTest = test_fields,
                                       defFieldTol = 1.e-4,
                                       expected = AUTOUGH2_result[run_name],
                                       testOutputIndex = -1))

jrunner = SimpleJobRunner(mpi = True)
testResult, mResults = co2_column_test.runTest(jrunner, createReports = True)

z = [lay.centre for lay in geo.layerlist[1:]]
for run_index, run_name in enumerate(run_names):
    for field_name in plot_fields:
        scale = field_scale[field_name]
        unit = field_unit[field_name]
        result = co2_column_test.mSuite.resultsList[run_index]
        var = result.getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(var, z, '-', label = 'Waiwera')
        var = AUTOUGH2_result[run_name].getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(var, z, 'o', label = 'AUTOUGH2')
        plt.xlabel(field_name + ' (' + unit + ')')
        plt.ylabel('elevation (m)')
        plt.title(' '.join(['Input CO2 mass fraction:', run_name + '%']))
        img_filename_base = '_'.join((model_name, run_name, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(co2_column_test.mSuite.runs[run_index].basePath,
                                    co2_column_test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.legend(loc = 'best')
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename)
        plt.clf()
        co2_column_test.mSuite.analysisImages.append(img_filename)

# generate report:

for rGen in getGenerators(["RST"], co2_column_test.outputPathBase):
    report_filename = os.path.join(co2_column_test.outputPathBase,
                     "%s-report.%s" % (co2_column_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(co2_column_test, mResults, rGen, report_filename)
    html_filename = os.path.join(co2_column_test.outputPathBase,
                     "%s-report.%s" % (co2_column_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")
