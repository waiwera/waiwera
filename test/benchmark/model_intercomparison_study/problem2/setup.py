from t2data import *
from t2incons import *
from t2thermo import cowat
from builtins import zip
from math import ceil
import json
import os

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'problem2'

# This is not the same mesh used in the original MIS problem, which
# specified cell centres that could not possibly be at the centres of
# the cells, and had a comparatively large cell next to the well. The
# mesh used here is chosen so that the second and third cell centres
# do correspond to cell centres in the original mesh.
dr = np.concatenate((np.array([0.3, 0.4, 0.6]),
                     0.6 * np.logspace(1, 29, 30, base = 1.28)))
thickness = 100.
dz = np.array([thickness])

geo = mulgrid().rectangular(dr, [1.], dz, atmos_type = 2)
geo.write('g' + model_name + '.dat')

dat = t2data()
dat.title = 'Model intercomparison study problem 2'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().radial(dr, dz)
rock = dat.grid.rocktypelist[0]
rock.porosity = 0.2
rock.permeability = np.ones(3) * 1.e-14
rock.density = 2650. # assumed
rock.conductivity = 0.0
rock.specific_heat = 2650. / rock.density * 1.e3

dat.relative_permeability = {'type': 1, 'parameters': [0., 0., 0., 0., 0.]}
dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}

ndt = 23
dts = 5. * np.logspace(0, ndt - 1, ndt, base = 1.5)
day = 24. * 60. * 60.

P0, T0 = 90.e5, 260.
Tinj = 160.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1. * day,
     'print_interval': 1,
     'gravity': None,
     'default_incons': [P0, T0],
     'const_timestep': -int(ceil(ndt / 8.)),
     'timestep': list(dts),
     'relative_error': 1.e-6,
     'absolute_error': 1.
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 1 # generation table interpolation
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 0 # initial output

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}

lay = geo.layerlist[-1]
col = geo.columnlist[0]
genrate = -14.
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = blkname, block = blkname,
                  gx = genrate, type = 'MASS')
dat.add_generator(gen)

dat.write(model_name + 'a.dat')

# incons:
inc = t2incon()
for blk in dat.grid.blocklist:
    inc[blk.name] = [P0, T0]
inc.write(model_name + 'a.incon')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + 'a.incon',
        silent = True)

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')
jsondata = dat.json(geo, mesh_filename, incons = inc, bdy_incons = inc,
                    mesh_coords = 'rz')
jsondata['initial']['primary'] = [P0, T0]
jsondata['output']['initial'] = False
jsondata['mesh']['radial'] = True
jsondata['output']['fields'] = {'fluid': ['liquid_saturation']}
json.dump(jsondata, file(model_name + 'a.json', 'w'), indent = 2)

# problem 2b:

P0 = 30.e5
S0 = 1. - 0.65
dat.parameter['default_incons'] = [P0, S0]
rock.porosity = 0.15
rock.permeability = np.ones(3) * 0.24e-12
rock.specific_heat = 2000. / rock.density * 1.e3
gen.gx = -16.7

dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.05, 0., 0., 0.]}

dat.write(model_name + 'b.dat')

inc = t2incon()
for blk in dat.grid.blocklist:
    inc[blk.name] = [P0, S0]
inc.write(model_name + 'b.incon')

dat.run(simulator = 'AUTOUGH2_41Da',
        incon_filename = model_name + 'b.incon',
        silent = True)

jsondata = dat.json(geo, mesh_filename, incons = inc, bdy_incons = inc,
                    mesh_coords = 'rz')
jsondata['initial'] = {'primary': [P0, S0], 'region': 4}
jsondata['output']['initial'] = False
jsondata['output']['fields'] = {'fluid': ['liquid_saturation']}
jsondata['mesh']['radial'] = True
json.dump(jsondata, file(model_name + 'b.json', 'w'), indent = 2)

# problem 2c:

P0 = 90.e5
T0 = 300.
dat.parameter['default_incons'] = [P0, T0]
rock.porosity = 0.2
rock.permeability = np.ones(3) * 0.01e-12
rock.specific_heat = 2650. / rock.density * 1.e3
gen.gx = -14.

dat.write(model_name + 'c.dat')

inc = t2incon()
for blk in dat.grid.blocklist:
    inc[blk.name] = [P0, T0]
inc.write(model_name + 'c.incon')

dat.run(simulator = 'AUTOUGH2_41Da',
        incon_filename = model_name + 'c.incon',
        silent = True)

jsondata = dat.json(geo, mesh_filename, incons = inc, bdy_incons = inc,
                    mesh_coords = 'rz')
jsondata['initial'] = {'primary': [P0, T0], 'region': 1}
jsondata['output']['initial'] = False
jsondata['output']['fields'] = {'fluid': ['liquid_saturation']}
jsondata['mesh']['radial'] = True
json.dump(jsondata, file(model_name + 'c.json', 'w'), indent = 2)

os.chdir(orig_dir)
