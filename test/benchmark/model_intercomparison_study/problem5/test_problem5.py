"""
Model Intercomparison Study problem 5
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

from credo.systest import SciBenchmarkTest
from credo.systest import FieldWithinTolTC
from credo.systest import HistoryWithinTolTC

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import ModelResult, HistoryDataResult
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from mulgrids import mulgrid
import matplotlib.pyplot as plt
import numpy as np

from docutils.core import publish_file

parser = argparse.ArgumentParser()
parser.add_argument("-np", type = int, default = 1, help = "number of processes")
parser.add_argument("-d", "--docker", action = "store_true",
                    help = "run via Docker (waiwera-dkr)")
args = parser.parse_args()
mpi = args.np > 1 and not args.docker
simulator = 'waiwera-dkr -np %d' % args.np if args.docker else 'waiwera'

def total_steam_history(mResult, index):
    """Returns history of total steam in place (per unit reservoir
    thickness) in the model."""
    thickness = 100.
    times = mResult.getTimes()
    steam_history = []
    for i in range(len(times)):
        sv = mResult.getFieldAtOutputIndex('Vapour saturation', i)
        rhov = mResult.getFieldAtOutputIndex('Vapour density', i)
        phi = mResult.getFieldAtOutputIndex('rock_porosity', i)
        vol = mResult.getFieldAtOutputIndex('geom_volume', i)
        steam = 1.e-3 * np.sum(phi * sv * rhov * vol) / thickness
        steam_history.append(steam)
    return times, np.array(steam_history)

model_name = 'problem5'

AUT2_FIELDMAP = {
    'Steam': total_steam_history
}
WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Vapour saturation': 'fluid_vapour_saturation',
    'Vapour density': 'fluid_vapour_density',
    'Enthalpy': 'source_enthalpy',
    'Steam': total_steam_history
}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

run_names = ['a', 'b']
test_fields = ["Pressure", "Temperature", "Vapour saturation"]
plot_fields = test_fields
digitised_test_fields = {'production': ["Pressure", "Temperature"],
                         'injection': ["Pressure"],
                         'total': ["Steam"]}
digitised_simulators = ["LBL", "S-Cubed"] 
test_source_fields = ["Enthalpy"]

geo = mulgrid(t2geo_filename)
map_out_atm = list(range(geo.num_atmosphere_blocks, geo.num_blocks))

problem5_test = SciBenchmarkTest(model_name + "_test", nproc = args.np)
problem5_test.description = """Model Intercomparison Study problem 5"""

def geo_pos_block_index(pos):
    # returns index of bottom layer block containing (x,y) pos
    col = geo.column_containing_point(pos)
    blkname = geo.block_name(geo.layerlist[-1].name, col.name)
    return geo.block_name_index[blkname]

obs_points = ["production", "injection"]
obs_positions = {"production": np.array([62.5, 62.5]),
             "injection": np.array([162.5, 137.5])}
obs_cell_index = dict([(obspt, geo_pos_block_index(obs_positions[obspt]))
                       for obspt in obs_points])
obs_points.append("total")
obs_cell_index["total"] = 0 # dummy value- no location

for run_index, run_name in enumerate(run_names):

    run_base_name = model_name + run_name
    run_filename = run_base_name + '.json'
    model_run = WaiweraModelRun(run_name, run_filename,
                              fieldname_map = WAIWERA_FIELDMAP,
                              simulator = simulator,
                              basePath = os.path.realpath(model_dir))
    model_run.jobParams['nproc'] = args.np
    problem5_test.mSuite.addRun(model_run, run_name)

problem5_test.setupEmptyTestCompsList()
digitised_result = {}
digitised_source_result = {}
AUTOUGH2_result = {}
source_index = 0

for run_index, run_name in enumerate(run_names):

    run_base_name = model_name + run_name
    dat_filename = os.path.join(model_dir, run_base_name + ".dat")
    run_filename = os.path.join(model_dir, run_base_name + ".listing")
    AUTOUGH2_result[run_name] = T2ModelResult("AUTOUGH2", run_filename,
                                              geo_filename = t2geo_filename,
                                              dat_filename = dat_filename,
                                              ordering_map = map_out_atm,
                                              fieldname_map = AUT2_FIELDMAP)

    for sim in digitised_simulators:
        data = {}
        for obspt in obs_points:
            for field_name in digitised_test_fields[obspt]:
                data_filename = '_'.join((model_name + run_name, obspt, field_name, sim))
                data_filename = os.path.join(data_dir, data_filename.lower() + '.dat')
                data[field_name, obs_cell_index[obspt]] = np.loadtxt(data_filename)
        digitised_result[run_name, sim] = HistoryDataResult(sim, data)
        obspt = "production"
        for field_name in test_source_fields:
            data_filename = '_'.join((model_name + run_name, obspt, field_name, sim))
            data_filename = os.path.join(data_dir, data_filename.lower() + '.dat')
            data[field_name, source_index] = np.loadtxt(data_filename)
        digitised_source_result[run_name, sim] = HistoryDataResult(sim, data)
    
    for obspt in obs_points:
        blk_index = obs_cell_index[obspt]
        if obspt == 'total': fields = digitised_test_fields[obspt]
        else: fields = test_fields
        problem5_test.addTestComp(run_index, "AUTOUGH2 " + obspt,
                              HistoryWithinTolTC(fieldsToTest = fields,
                                                 defFieldTol = 1.e-3,
                                                 expected = AUTOUGH2_result[run_name],
                                                 testCellIndex = blk_index))
        for sim in digitised_simulators:
            problem5_test.addTestComp(run_index, ' '.join((sim, obspt)),
                                      HistoryWithinTolTC(fieldsToTest = \
                                                         digitised_test_fields[obspt],
                                                         defFieldTol = 1.5e-2,
                                                         fieldTols = {"Steam": 5.e-2},
                                                         expected = digitised_result[run_name, sim],
                                                         testCellIndex = obs_cell_index[obspt],
                                                         orthogonalError = True))
    obspt = "production"
    problem5_test.addTestComp(run_index, "AUTOUGH2 " + obspt + " source",
                              HistoryWithinTolTC(fieldsToTest = test_source_fields,
                                                 defFieldTol = 1.e-3,
                                                 expected = AUTOUGH2_result[run_name],
                                                 testSourceIndex = source_index))
    for sim in digitised_simulators:
        problem5_test.addTestComp(run_index, ' '.join((sim, obspt, "source")),
                                  HistoryWithinTolTC(fieldsToTest = \
                                                     test_source_fields,
                                                     defFieldTol = 1.5e-2,
                                                     expected = digitised_source_result[run_name, sim],
                                                     testSourceIndex = source_index,
                                                     orthogonalError = True))

jrunner = SimpleJobRunner(mpi = mpi)
testResult, mResults = problem5_test.runTest(jrunner, createReports = True)

# plots:
elevations = [layer.centre for layer in geo.layerlist[1:]]
scale = {"Pressure": 1.e5, "Temperature": 1., "Steam": 1., "Enthalpy": 1.e3}
unit = {"Pressure": "bar", "Temperature": "$^{\circ}$C", "Steam": "tonnes/m",
        "Enthalpy": "kJ/kg"}
symbol = {"LBL": 's', "S-Cubed": 'o'}
yr = 365. * 24. * 60. * 60.

for run_index, run_name in enumerate(run_names):
    for obspt in obs_points:

        # plot time history results:
        tc_name = "AUTOUGH2 " + obspt
        blk_index = obs_cell_index[obspt]

        for field_name in digitised_test_fields[obspt]:

            t, var = problem5_test.mSuite.resultsList[run_index].\
                     getFieldHistoryAtCell(field_name, blk_index)
            plt.plot(t / yr, var / scale[field_name], '-', label = 'Waiwera')

            t, var = AUTOUGH2_result[run_name].getFieldHistoryAtCell(field_name, blk_index)
            plt.plot(t / yr, var / scale[field_name], '+', label = 'AUTOUGH2')

            for sim in digitised_simulators:
                result = digitised_result[run_name, sim]
                t, var = result.getFieldHistoryAtCell(field_name, blk_index)
                plt.plot(t / yr, var / scale[field_name], symbol[sim], label = sim)
            plt.xlabel('time (years)')
            plt.ylabel(field_name + ' (' + unit[field_name] + ')')
            plt.legend()
            plt.title(' '.join((model_name + run_name, field_name.lower(),
                                'results:', obspt)))
            img_filename_base = '_'.join((model_name, run_name, obspt, field_name))
            img_filename_base = img_filename_base.replace(' ', '_')
            img_filename = os.path.join(problem5_test.mSuite.runs[run_index].basePath,
                                        problem5_test.mSuite.outputPathBase,
                                        img_filename_base + '.png')
            plt.tight_layout(pad = 3.)
            plt.savefig(img_filename)
            plt.clf()
            problem5_test.mSuite.analysisImages.append(img_filename)

            if obspt != 'total':
                t = problem5_test.testComps[run_index][tc_name].times
                for field_name in plot_fields:
                    var = np.array(problem5_test.testComps[run_index][tc_name].fieldErrors[field_name])
                    plt.plot(t / yr, var, '-o')
                    plt.xlabel('time (years)')
                    plt.ylabel(field_name + ' error')
                    plt.title(' '.join((model_name + run_name, 'comparison with AUTOUGH2: ',
                                       obspt)))
                    img_filename_base = '_'.join((model_name, run_name, tc_name, 'error', field_name))
                    img_filename_base = img_filename_base.replace(' ', '_')
                    img_filename = os.path.join(problem5_test.mSuite.runs[run_index].basePath,
                                                problem5_test.mSuite.outputPathBase,
                                                img_filename_base + '.png')
                    plt.tight_layout(pad = 3.)
                    plt.savefig(img_filename)
                    plt.clf()
                    problem5_test.mSuite.analysisImages.append(img_filename)

    obspt = "production"
    source_index = 0
    for field_name in test_source_fields:

        t, var = problem5_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t / yr, var / scale[field_name], '-', label = 'Waiwera')

        t, var = AUTOUGH2_result[run_name].getFieldHistoryAtSource(field_name, source_index)
        plt.plot(t / yr, var / scale[field_name], '+', label = 'AUTOUGH2')

        for sim in digitised_simulators:
            result = digitised_source_result[run_name, sim]
            t, var = result.getFieldHistoryAtSource(field_name, source_index)
            plt.plot(t / yr, var / scale[field_name], symbol[sim], label = sim)
        plt.xlabel('time (years)')
        plt.ylabel(field_name + ' (' + unit[field_name] + ')')
        plt.legend()
        plt.title(' '.join((model_name + run_name, field_name.lower(),
                            'results:', obspt)))
        img_filename_base = '_'.join((model_name, run_name, obspt, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(problem5_test.mSuite.runs[run_index].basePath,
                                    problem5_test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename)
        plt.clf()
        problem5_test.mSuite.analysisImages.append(img_filename)

# generate report:

for rGen in getGenerators(["RST"], problem5_test.outputPathBase):
    report_filename = os.path.join(problem5_test.outputPathBase,
                     "%s-report.%s" % (problem5_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(problem5_test, mResults, rGen, report_filename)
    html_filename = os.path.join(problem5_test.outputPathBase,
                     "%s-report.%s" % (problem5_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")

