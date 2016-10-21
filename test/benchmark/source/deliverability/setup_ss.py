# Set up steady state solutions for AUTOUGH2 and supermodel, to be
# used as initial conditions for transient runs.

import os
from t2data_json import *
from t2thermo import *
import json
import subprocess

model_name = 'deliv'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

t2geo_filename = 'g' + model_name + '.dat'
t2dat_filename = model_name + '_ss.dat'

d = 100.
dx = [d]
dy = [d]
dz = [d] * 10

geo = mulgrid().rectangular(dx, dy, dz, atmos_type = 0)
geo.write(t2geo_filename)

pressure = 1.e5
temperature = 20.

dat = t2data_export_json()
dat.title = 'Deliverability source test'
dat.simulator = 'AUTOUGH2.2'
dat.grid = t2grid().fromgeo(geo)
dat.grid.rocktype['dfalt'].permeability = np.ones(3) * 1.e-13

ndt = 200
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e15,
     'print_interval': ndt,
     'gravity': 9.8,
     'default_incons': [pressure, temperature],
     'const_timestep': 1.e4
     })
dat.start = True
dat.parameter['option'][1] = 1
dat.parameter['option'][16] = 5

dat.multi = {'num_components':1, 'num_equations':2,
             'num_phases':2, 'num_secondary_parameters':6,
             'eos': 'EW'}

dat.relative_permeability = {'type': 1, 'parameters': [0.,0.,1., 1.]}
dat.capillarity = {'type': 1, 'parameters': [0.,0.,0., 0.]}

gx = 10. # kg/s
t = 240. # deg C
p = sat(t)
d,u = cowat(t,p)
h = u + p / d

gen = t2generator(name = 'gen 1', block = dat.grid.blocklist[-1].name,
                  gx = gx, ex = h, type = 'MASS')
dat.add_generator(gen)
dat.write(t2dat_filename)

# run steady state
t2gener_filename = 'gener.data'
if os.path.isfile(t2gener_filename): os.remove(t2gener_filename)
dat.run(simulator = 'AUTOUGH2_41Dasw', silent = True)

mesh_filename = 'g'+ model_name +'.exo'
geo.write_exodusii(mesh_filename)
jsondata = dat.json(geo, mesh_filename)
jsondata["output"]["frequency"] = 0
jsondata["output"]["initial"] = False
filename = model_name + '_ss.json'
json.dump(jsondata, file(filename, 'w'), indent = 2)
subprocess.check_output(['mpiexec', '-np', '1', 'supermodel', filename])

os.chdir(orig_dir)
