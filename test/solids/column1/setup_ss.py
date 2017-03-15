# Hot injection into cold column- steady state

from t2data_json import *
import json

name = 'hotcol'
d = 100.
dx = [d]
dy = [d]
dz = [d] * 10

geo = mulgrid().rectangular(dx, dy, dz, atmos_type = 0)
geo.write('g'+ name +'.dat')

pressure, temperature = 1.e5, 20.

dat = t2data_export_json()
dat.title = 'Hot injection into cold column'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)
dat.grid.rocktype['dfalt'].permeability = np.ones(3) * 1.e-13

ndt = 200
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e14,
     'print_interval': ndt,
     'gravity': 9.8,
     'default_incons': [pressure, temperature],
     'const_timestep': 1.e4,
     'relative_error': 1.e-6,
     'absolute_error': 1.
     })
dat.start = True
dat.parameter['option'][1] = 1
dat.parameter['option'][16] = 6
dat.parameter['option'][21] = 1 # direct solver

dat.multi = {'num_components':1, 'num_equations':2,
             'num_phases':2, 'num_secondary_parameters':6,
             'eos': 'EW'}

dat.relative_permeability = {'type': 1, 'parameters': [0.,0.,1., 1.]}
dat.capillarity = {'type': 1, 'parameters': [0.,0.,0., 0.]}

# run steady state
from os import remove
from os.path import isfile
if isfile('gener.data'): remove('gener.data')
dat.write(name + '_ss.dat')
dat.run(simulator = 'AUTOUGH2_41D')

# export input to Waiwera
mesh_filename = 'g'+ name +'.exo'
geo.write_exodusii(mesh_filename)
jsondata = dat.json(geo, mesh_filename)
jsondata["output"]["frequency"] = 0
jsondata["output"]["initial"] = False
filename = name + '_ss.json'
json.dump(jsondata, file(filename, 'w'), indent = 2)
