"""
CO2 column test from O'Sullivan et al. (1985)
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

from credo.systest import FieldWithinTolTC, OneDSolutionWithinTolTC

from mulgrids import mulgrid
from numpy import polyval

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
    'CO2 mass fraction': 'CO2 mass fractio',
    'CO2 partial pressure': 'CO2 partial pres'}

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'CO2 mass fraction': total_CO2_mass_fraction,
    'CO2 partial pressure': 'fluid_CO2_partial_pressure'}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')
geo = mulgrid(t2geo_filename)
map_out_atm = list(range(geo.num_atmosphere_blocks, geo.num_blocks))

CO2_mass_fractions = [0, 0.1, 1, 5] # percent
run_names = [str(xgp) for xgp in CO2_mass_fractions]
obs_cell_index = 0

test_fields = ['Pressure', 'Temperature', 'Vapour saturation', 'CO2 mass fraction']
plot_fields = ['Vapour saturation', 'CO2 mass fraction', 'CO2 partial pressure']
field_scale = {'Temperature': 1., 'Vapour saturation': 1., 'CO2 mass fraction': 1.,
               'CO2 partial pressure': 1.e5}
field_unit = {'Temperature': '$^{\circ}$C', 'Vapour saturation': '', 'CO2 mass fraction': '',
              'CO2 partial pressure': 'bar'}

digitised_test_fields = {'CO2 partial pressure'}
digitised_simulators = ['MULKOM']
digitised_run_names = ['0.1', '1', '5']

co2_column_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
co2_column_test.description = """Vertical column CO2 test from O'Sullivan et al. (1985). The original
problem (figs 10 and 11) specified CO2 input in terms of partial pressure, though the text indicated
that mass fraction (%) may actually have been used. Here, CO2 input is in terms of mass fraction. The
vapour saturation results, even for AUTOUGH2, do not match those in the original paper; however,
the problem is still useful as a benchmark of Waiwera against AUTOUGH2.
"""

for run_index, run_name in enumerate(run_names):
    run_base_name = model_name + '_' + run_name
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                                fieldname_map = WAIWERA_FIELDMAP,
                                simulator = simulator,
                                basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = args.np
    co2_column_test.mSuite.addRun(model_run, run_name)

co2_column_test.setupEmptyTestCompsList()

digitised_test_fields = ['CO2 partial pressure']
digitised_simulators = ['MULKOM']
digitised_result = {}
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
                                       defFieldTol = 1.e-3,
                                       expected = AUTOUGH2_result[run_name],
                                       testOutputIndex = -1))

for run_index, run_name in enumerate(run_names):
    if run_name in digitised_run_names:
        for sim in digitised_simulators:
            for field_name in digitised_test_fields:
                data_filename = '_'.join((model_name, run_name, field_name, sim))
                data_filename = data_filename.lower().replace(' ', '_')
                data_filename = os.path.join(data_dir, data_filename + '.dat')
                result = DigitisedOneDFieldResult(sim, data_filename, field_name, -1)
                result.data[:,1] *= field_scale[field_name]
                digitised_result[run_name, field_name, sim] = result
                co2_column_test.addTestComp(run_index, ' '.join((sim, field_name)),
                                          OneDSolutionWithinTolTC(
                                              fieldsToTest = [field_name],
                                              defFieldTol = 0.07,
                                              expected = result,
                                              coordinateIndex = 1,
                                              testOutputIndex = -1))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = co2_column_test.runTest(jrunner, createReports = True)

z = [lay.centre for lay in geo.layerlist[1:]]
for run_index, run_name in enumerate(run_names):
    for field_name in plot_fields:
        scale = field_scale[field_name]
        unit = field_unit[field_name]
        result = co2_column_test.mSuite.resultsList[run_index]
        var = result.getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(var, z, '-', label = 'Waiwera', zorder = 2)
        var = AUTOUGH2_result[run_name].getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(var, z, 's', label = 'AUTOUGH2', zorder = 1)
        if field_name in digitised_test_fields and run_name in digitised_run_names:
            for sim in digitised_simulators:
                result = digitised_result[run_name, field_name, sim]
                zd = result.getCoordinates()
                var = result.getFieldAtOutputIndex(field_name, -1) / scale
                plt.plot(var, zd, 'o', label = sim)
        plt.xlabel(field_name + ' (' + unit + ')')
        plt.ylabel('elevation (m)')
        plt.title(' '.join(['Input CO2 mass fraction:', run_name + '%']))
        img_filename_base = '_'.join((model_name, run_name, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(co2_column_test.mSuite.runs[run_index].basePath,
                                    co2_column_test.mSuite.outputPathBase,
                                    img_filename_base)
        plt.legend(loc = 'best')
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename + '.png', dpi = 300)
        plt.savefig(img_filename + '.pdf')
        plt.clf()
        co2_column_test.mSuite.analysisImages.append(img_filename + '.png')

# combined partial pressure plot:
field_name = 'CO2 partial pressure'
for run_index, run_name in enumerate(run_names):
    if run_name in digitised_run_names:
        scale = field_scale[field_name]
        unit = field_unit[field_name]
        result = co2_column_test.mSuite.resultsList[run_index]
        var = result.getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(var, z, 'b-', label = 'Waiwera', zorder = 2)
        var = AUTOUGH2_result[run_name].getFieldAtOutputIndex(field_name, -1) / scale
        plt.plot(var, z, 'gs', label = 'AUTOUGH2', zorder = 1)
        if field_name in digitised_test_fields:
            for sim in digitised_simulators:
                result = digitised_result[run_name, field_name, sim]
                zd = result.getCoordinates()
                var = result.getFieldAtOutputIndex(field_name, -1) / scale
                plt.plot(var, zd, 'ro', label = sim)
plt.text(2.3, -900, "0.1%")
plt.text(10, -800, "1%")
plt.text(18, -800, "5%")
plt.xlabel(field_name + ' (' + unit + ')')
plt.ylabel('elevation (m)')
img_filename_base = '_'.join((model_name, field_name))
img_filename_base = img_filename_base.replace(' ', '_')
img_filename = os.path.join(co2_column_test.mSuite.runs[run_index].basePath,
                            co2_column_test.mSuite.outputPathBase,
                            img_filename_base)
# remove duplicate labels in legend:
handles, labels = plt.gca().get_legend_handles_labels()
i = 1
while i < len(labels):
    if labels[i] in labels[:i]:
        del(labels[i])
        del(handles[i])
    else: i += 1
plt.legend(handles, labels, loc = 'best')
plt.tight_layout(pad = 3.)
plt.savefig(img_filename + '.png', dpi = 300)
plt.savefig(img_filename + '.pdf')
plt.clf()

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
