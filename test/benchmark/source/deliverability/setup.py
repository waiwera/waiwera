# Transient model run setup.
# Before running this script, setup_ss.py must be run to compute the steady state solutions.

import os
from t2data_json import *
from t2thermo import *
from t2incons import *
import json

model_name = 'deliv'

model_dir = './run'
orig_dir = os.getcwd()
os.chdir(model_dir)

inc = t2incon(model_name + '_ss.save')
inc.write(model_name + '.incon')

geo = mulgrid('g' + model_name + '.dat')
dat = t2data_export_json(model_name + '_ss.dat')

mesh_filename = 'g'+ model_name +'.exo'

ndt = 150
maxtime = 1.e9
PI = 1.e-11
Pwb = 5.e5

dat.start = False

dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': maxtime,
     'print_interval': 1,
     'relative_error': 1.e-7,
     })

dat.parameter['option'][24] = 2 # output initial results
dat.parameter['option'][16] = 0
dat.parameter['timestep'] = list(np.logspace(4, 7.5, 80))
dat.parameter['const_timestep'] = -int(len(dat.parameter['timestep']) / 8 + 0.5)

dat.write(model_name + '.dat')

iblk = geo.num_blocks / 2

run_names = ['delv', 'delg_flow', 'delg_pi_table', 'delg_limit', 'delt', 'delw']
run_sources = {}

# straight DELV
run_sources['delv'] = t2generator(type = 'DELV', block = geo.block_name_list[iblk],
                                  name = 'del 1', gx = PI, ex = Pwb)

# DELG with PI calculated from flow rate
run_sources['delg_flow'] = t2generator(type = 'DELG', block = geo.block_name_list[iblk],
                                       name = 'del 1', ex = Pwb, hg = -20.)

# DELG with PI table vs time
run_sources['delg_pi_table'] = t2generator(type = 'DELG', block = geo.block_name_list[iblk],
                                           name = 'del 1', ex = Pwb, ltab = 4,
                                           time = [0., 1.e4, 1.e6, 1.e9],
                                           rate = [1.e-11, 0.7e-11, 0.2e-11, 0.1e-11])

# DELG with steam limiter
run_sources['delg_limit'] = t2generator(type = 'DELG', block = geo.block_name_list[iblk],
                                        name = 'del 1', gx = PI, ex = Pwb, hg = 2.5, fg = 0.7e6)

# DELT with total limiter
run_sources['delt'] = t2generator(type = 'DELT', block = geo.block_name_list[iblk],
                                  name = 'del 1', gx = PI, ex = Pwb, hg = 20., fg = 0.7e6)

# DELW with water limiter
run_sources['delw'] = t2generator(type = 'DELW', block = geo.block_name_list[iblk],
                     name = 'del 1', gx = PI, ex = Pwb, hg = 20., fg = 0.7e6)

from os import remove
from os.path import isfile

for run_name in run_names:

    dat = t2data_export_json(model_name + '.dat')
    dat.add_generator(run_sources[run_name])
    run_base_name = model_name + '_' + run_name
    dat.write(run_base_name + '.dat')

    if isfile('gener.data'): remove('gener.data')
    dat.run(simulator = 'AUTOUGH2_41Dasw',
            incon_filename = model_name + '.incon',
            silent = True)

    jsondata = dat.json(geo, mesh_filename, incons = model_name + '_ss.h5',
                    bdy_incons = inc)
    filename = run_base_name + '.json'
    json.dump(jsondata, file(filename, 'w'), indent = 2)

os.chdir(orig_dir)
