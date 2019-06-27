# One-cell CO2 test from O'Sullivan et al. (1985)
from __future__ import print_function
from t2data import *
from t2incons import *
from t2thermo import *
import json
import os

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'co2_one_cell'

dx = 1.
geo = mulgrid().rectangular([dx], [dx], [dx], atmos_type = 2)
geo.write('g' + model_name + '.dat')

dat = t2data()
dat.title = 'CO2 one cell'
dat.simulator = 'AUTOUGH2.2EWC'

dat.grid = t2grid().fromgeo(geo)

rock = dat.grid.rocktypelist[0]
rock.density = 2500.
rock.porosity = 0.15
rock.specific_heat = 900.

# The exact Corey parameters are not stated in the paper and have to
# be estimated from fig 4:
dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.05, 0., 0., 0.]}
dat.capillarity = {'type': 8, 'parameters': [0., 0., 0., 0., 0.]}

dat.multi = {'eos': 'EWC', 'num_components': 2, 'num_phases': 2,
             'num_equations': 3, 'num_secondary_parameters': 6}
dat.start = True

ndt = 50
T0, Sg0, PCO20 = 260., 0.2, 30.e5
P0 = sat(T0) + PCO20
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 19.,
     'print_interval': 1,
     'gravity': None,
     'default_incons': [P0, Sg0, PCO20],
     'const_timestep': 0.5,
     'max_timestep': 0.5
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 2 # initial output

gen = t2generator(block = dat.grid.blocklist[0].name,
                  name = 'gen 1', type = 'MASS', gx = -5.)
dat.add_generator(gen)

dat.write(model_name + '.dat')
dat.run(simulator = AUTOUGH2, silent = True)

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')
jsondata = dat.json(geo, mesh_filename, mesh_coords = 'xz')
jsondata['initial']['primary'] = [P0, Sg0, PCO20]
jsondata['initial']['region'] = 4
jsondata['gravity'] = None
json.dump(jsondata, file(model_name + '.json', 'w'), indent = 2)

os.chdir(orig_dir)
