# Set up deliverability benchmark test cases

import os
from t2data_json import *
from t2thermo import *
import json

AUTOUGH2 = 'AUTOUGH2_42D'

model_name = 'deliv'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

t2geo_filename = 'g' + model_name + '.dat'
t2dat_filename = model_name + '.dat'

d = 100.
dx = [d] * 10
dy = [d]
dz = [d]

geo = mulgrid().rectangular(dx, dy, dz, atmos_type = 2)
geo.write(t2geo_filename)

dp = 1.e5
temperature = 230.
pressure = sat(temperature) + dp

dat = t2data_export_json()
dat.title = 'Deliverability source test'
dat.simulator = 'AUTOUGH2.2'
dat.grid = t2grid().fromgeo(geo)
dat.grid.rocktype['dfalt'].permeability = np.ones(3) * 1.e-13

ndt = 80
dts = np.logspace(4, 7.5, ndt)
maxtime = np.sum(dts)

dat.start = False

dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': maxtime,
     'print_interval': 1,
     'default_incons': [pressure, temperature],
     'relative_error': 1.e-7,
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][16] = 0

dat.parameter['timestep'] = list(dts)
dat.parameter['const_timestep'] = -int(round(ndt // 8))

dat.start = True

dat.multi = {'num_components':1, 'num_equations':2,
             'num_phases':2, 'num_secondary_parameters':6,
             'eos': 'EW'}

dat.relative_permeability = {'type': 1, 'parameters': [0.,0.,1., 1.]}
dat.capillarity = {'type': 1, 'parameters': [0.,0.,0., 0.]}

# BC block at right hand end:
x = 1.e3
col = geo.columnlist[-1]
lay = geo.layerlist[-1]
blkname = geo.block_name(lay.name, col.name)
blk = dat.grid.block[blkname]
vol = 0.
rock = dat.grid.rocktypelist[-1]
bcblk = t2block('bc %2d' % 0, vol, rock,
                np.array([x, 0.5 * d, -0.5 * d]))
dat.grid.add_block(bcblk)
area = d * d
con = t2connection([blk, bcblk], 1,
                   [0.5 * d, geo.atmosphere_connection],
                   area, 0.)
dat.grid.add_connection(con)

dat.write(t2dat_filename)

run_names = ['delv', 'delg_flow', 'delg_pi_table', 'delg_pwb_table', 'delg_limit', 'delt', 'delw']
run_sources = {}
iblk = 0
PI = 1.e-11
Pwb = 5.e5

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

# DELG with Pwb table vs enthalpy
run_sources['delg_pwb_table'] = t2generator(type = 'DELG', block = geo.block_name_list[iblk],
                                            name = 'del 1', gx = PI, ltab = -3,
                                            time = [0., 1.1e6, 2.8e6],
                                            rate = [22.e5, 20.e5, 0.5e5])

# DELG with steam limiter
run_sources['delg_limit'] = t2generator(type = 'DELG', block = geo.block_name_list[iblk],
                                        name = 'del 1', gx = PI, ex = Pwb, hg = 2.5, fg = 0.7e6)

# DELT with total limiter
run_sources['delt'] = t2generator(type = 'DELT', block = geo.block_name_list[iblk],
                                  name = 'del 1', gx = PI, ex = Pwb, hg = 20., fg = 0.7e6)

# DELW with water limiter
run_sources['delw'] = t2generator(type = 'DELW', block = geo.block_name_list[iblk],
                     name = 'del 1', gx = PI, ex = Pwb, hg = 20., fg = 0.7e6)

mesh_filename = 'g' + model_name + '.exo'
geo.write_mesh(mesh_filename)

from os import remove
from os.path import isfile

for run_name in run_names:

    dat = t2data_export_json(model_name + '.dat')
    if run_name == 'delg_pi_table':
        # AUTOUGH2 initial flow rate appears to be wrong for this case:
        dat.parameter['option'][24] = 0
    else: # output initial results
        dat.parameter['option'][24] = 2
    dat.add_generator(run_sources[run_name])
    run_base_name = model_name + '_' + run_name
    dat.write(run_base_name + '.dat')

    if isfile('gener.data'): remove('gener.data')
    dat.run(simulator = AUTOUGH2,
            incon_filename = model_name + '.incon',
            silent = True)

    jsondata = dat.json(geo, mesh_filename)
    filename = run_base_name + '.json'
    json.dump(jsondata, file(filename, 'w'), indent = 2)

os.chdir(orig_dir)
