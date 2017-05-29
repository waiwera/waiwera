# Set up AUTOUGH2 and Waiwera runs

import os
from t2data_json import *
from t2thermo import *
import json
import subprocess

AUTOUGH2 = 'AUTOUGH2_42D'

model_name = 'recharge'
run_name = 'outflow'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

t2geo_filename = 'g' + model_name + '.dat'
t2dat_filename = model_name + '_' + run_name + '.dat'

d = 100.
dx = [d] * 10
dy = [d]
dz = [d]

geo = mulgrid().rectangular(dx, dy, dz, atmos_type = 2)
geo.write(t2geo_filename)

pressure = 2.e5
temperature = 20.

dat = t2data_export_json()
dat.title = 'Recharge source test'
dat.simulator = 'AUTOUGH2.2'
dat.grid = t2grid().fromgeo(geo)
dat.grid.rocktype['dfalt'].permeability = np.ones(3) * 1.e-13

ndt = 25
dts = np.logspace(4, 6, ndt)
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': np.sum(dts),
     'print_interval': 1,
     'default_incons': [pressure, temperature],
     'relative_error': 1.e-7
     })
dat.start = True
dat.parameter['option'][1] = 1
dat.parameter['option'][24] = 2 # output initial results
dat.parameter['option'][16] = 0
dat.parameter['timestep'] = list(dts)
dat.parameter['const_timestep'] = -int(round(ndt // 8))

dat.multi = {'num_components':1, 'num_equations':2,
             'num_phases':2, 'num_secondary_parameters':6,
             'eos': 'EW'}

dat.relative_permeability = {'type': 1, 'parameters': [0.,0.,1., 1.]}
dat.capillarity = {'type': 1, 'parameters': [0.,0.,0., 0.]}

rech = t2generator(name = 'rec 1', block = dat.grid.blocklist[0].name,
                   gx = 1.e-3, fg = -1., hg = 1.0e5, type = 'RECH')
dat.add_generator(rech)

dat.write(t2dat_filename)

t2gener_filename = 'gener.data'
if os.path.isfile(t2gener_filename): os.remove(t2gener_filename)
dat.run(simulator = AUTOUGH2, silent = True)

mesh_filename = 'g'+ model_name +'.exo'
geo.write_exodusii(mesh_filename)
jsondata = dat.json(geo, mesh_filename)
jsondata["output"]["frequency"] = 1
jsondata["output"]["initial"] = True
filename = model_name + '_' + run_name + '.json'
json.dump(jsondata, file(filename, 'w'), indent = 2)

os.chdir(orig_dir)
