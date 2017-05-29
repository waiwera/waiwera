from t2data_json import *
from t2incons import *
from t2thermo import cowat
from builtins import zip
import json
import os

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'problem1'
dimensions = [1000., 100., 100.]
b = 0.5 * dimensions[1]
nblks = [40, 1, 1]
gridsizes = [dim // nblk for dim, nblk in zip(dimensions, nblks)]
dx = [[gridsize] * nblk for gridsize, nblk in zip(gridsizes, nblks)]

geo = mulgrid().rectangular(dx[0], dx[1], dx[2], atmos_type = 2)
geo.write('g' + model_name + '.dat')

dat = t2data_export_json()
dat.title = 'Model intercomparison study problem 1'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().radial(dx[0], dx[2])
rock = dat.grid.rocktypelist[0]
rock.porosity = 0.2
rock.density = 2500.
rock.permeability = np.ones(3) * 1.e-12
rock.conductivity = 20.0
rock.specific_heat = 1000.

dat.relative_permeability = {'type': 1, 'parameters': [0.] * 5}
dat.capillarity = {'type': 8, 'parameters': [0.] * 5}

ndt = 500
P0, T0 = 50.e5, 170.
Tinj = 160.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e9,
     'print_interval': 1,
     'gravity': None,
     'default_incons': [P0, T0],
     'const_timestep': -2.0,
     'timestep': [
         1.0000e5, 1.5000e5, 2.2500e5, 3.3750e5, 5.0625e5,
         7.5938e5, 1.1391e6, 1.7086e6, 2.5629e6, 3.8443e6,
         5.7665e6, 8.6498e6, 1.2975e7, 1.6700e7, 1.6700e7,
         1.6700e7],
     'max_timestep': 1.67e7,
     'relative_error': 1.e-6,
     'absolute_error': 1.
     })


dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 1 # generation table interpolation
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 2 # initial output

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}

lay = geo.layerlist[-1]
col = geo.columnlist[0]
genrate = 10.0
d,u = cowat(Tinj, P0)
h = u + P0 / d
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = blkname, block = blkname,
                  gx = genrate, ex = h, type = 'MASS')
dat.add_generator(gen)

# BC block at r = 1000:
col = geo.columnlist[-1]
blkname = geo.block_name(lay.name, col.name)
blk = dat.grid.block[blkname]
bcblk = t2block('bc %2d' % 0, geo.atmosphere_volume, rock,
                np.array([dimensions[0], 0., -0.5 * dimensions[2]]))
dat.grid.add_block(bcblk)
area = 2. * np.pi * dimensions[0] * gridsizes[2]
con = t2connection([blk, bcblk], 1,
                   [0.5 * gridsizes[0], geo.atmosphere_connection],
                   area, 0.)
dat.grid.add_connection(con)
dat.write(model_name + '.dat')

# incons:
inc = t2incon()
for blk in dat.grid.blocklist:
    inc[blk.name] = [P0, T0]
inc.write(model_name + '.incon')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '.incon',
        silent = True)

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')
jsondata = dat.json(geo, mesh_filename, incons = inc, bdy_incons = inc,
                    mesh_coords = 'rz')
jsondata['initial']['primary'] = [P0, T0]
jsondata['mesh']['radial'] = True
jsondata['gravity'] = [0., -9.8]
json.dump(jsondata, file(model_name + '.json', 'w'), indent = 2)
                  
os.chdir(orig_dir)
