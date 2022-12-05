from t2data import *
from t2incons import *
from builtins import zip
import json
import os

AUTOUGH2 = 'AUTOUGH2'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'salt_production'

rw = 1.
dr = np.logspace(0, 3, 40)
thickness = 500.
dz = np.array([thickness])

geo = mulgrid().rectangular(dr, [1.], dz, atmos_type = 2, origin = [rw, 0., 0.])
geo.write('g' + model_name + '.dat')

dat = t2data()
dat.title = 'Salt production problem'
dat.simulator = 'AUTOUGH2.2EWSG'

dat.grid = t2grid().radial(dr, dz, origin = [rw, 0.])
rock = dat.grid.rocktypelist[0]
rock.porosity = 0.05
rock.permeability = np.ones(3) * 50.e-15
rock.density = 2600.
rock.conductivity = 2.5
rock.specific_heat = 1000.

dat.relative_permeability = {'type': 1, 'parameters': [0.35, 0., 1., 0.7, 0.]}
dat.capillarity = {'type': 8, 'parameters': [0., 0., 0., 0., 0.]}

P0 = 43.23775e5
T0 = 275.55
Xs0 = 0.35

dat.parameter.update(
    {'max_timesteps': 100,
     'tstop': 1.e6,
     'print_interval': 1,
     'gravity': None,
     'const_timestep': 1.e3,
     'relative_error': 1.e-6,
     'absolute_error': 1.
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 1 # generation table interpolation
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 2 # initial output

lay = geo.layerlist[-1]
col = geo.columnlist[0]
genrate = -65.
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = 'wel 1', block = blkname,
                  gx = genrate, type = 'MASS')
dat.add_generator(gen)

dat.multi = {'eos': 'EWSG', 'num_components': 3, 'num_phases': 3,
             'num_equations': 4, 'num_secondary_parameters': 6}

dat.selection['integer'] = [0] * 9 + [0] + [0] + [0] * 5 + [2]
dat.parameter.update({'default_incons': [P0, Xs0, 0., T0]})

dat.write(model_name + '.dat')

inc = t2incon()
for blk in dat.grid.blocklist:
    inc[blk.name] = [P0, Xs0, 0., T0]
inc.write(model_name + '.incon')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '.incon',
        silent = True)

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x', file_format = 'gmsh22')

dat.simulator = 'AUTOUGH2.2EW'
dat.multi['eos'] = 'EW'
jsondata = dat.json(geo, mesh_filename, incons = inc, mesh_coords = 'rz')
jsondata['output']['initial'] = True
jsondata['mesh']['radial'] = True
jsondata['output'] = {'fields': {'fluid': ['liquid_saturation']}}

jsondata['eos'] = 'wse'
jsondata['initial']['primary'] = [P0, T0, Xs0]
jsondata['initial']['region'] = 1
json.dump(jsondata, open(model_name + '.json', 'w'), indent = 2)

os.chdir(orig_dir)
