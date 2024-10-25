# Near-critical 1D horizontal problem (Kissling, 2004)
from __future__ import print_function
from t2data import *
from t2incons import *
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

AUTOUGH2 = 'AUTOUGH2_2S'
model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'horizontal'

dx, nx = 10, 100
dy, dz = 10, 10
geo = mulgrid().rectangular([dx] * nx, [dy], [dz], atmos_type = 2)
geo.write('g' + model_name + '.dat')

dat = t2data()
dat.title = 'Near-critical 1-D horizontal problem'
dat.simulator = 'AUTOUGH2.2EW'

dat.grid = t2grid().fromgeo(geo)

rock = dat.grid.rocktypelist[0]
rock.density = 2650.
rock.porosity = 0.01
rock.permeability = np.ones(3) * 5e-15
rock.conductivity = 3.
rock.specific_heat = 1000.

dat.relative_permeability = {'type': 1, 'parameters': [0.3, 0.3, 0.7, 0.7, 0.]}
dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}
dat.start = True
dat.output_times = {'num_times_specified': 4, 'time': [1e6, 1e7, 1e8, 1e9]}

ndt = 1000
P0, T0 = 22.5e6, 374.

dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e11,
     'print_interval': ndt,
     'gravity': None,
     'default_incons': [P0, T0],
     'const_timestep': 1e4,
     'max_timestep': 1e11,
     'relative_error': 1.e-5,
     'absolute_error': 1.,
     'derivative_increment': 1.e-8
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][16] = 5 
dat.parameter['option'][19] = 1 # alternative initial conditions
dat.parameter['option'][20] = 1 # VPL off
dat.parameter['option'][23] = 0 # MULKOM compatibility off
dat.parameter['option'][24] = 2 # initial output

# BC block at x = 0:
col = geo.columnlist[0]
lay = geo.layerlist[-1]
blkname = geo.block_name(lay.name, col.name)
blk = dat.grid.block[blkname]
vol = 0.
bc_name = 'bc  0'
bcblk = t2block(bc_name, vol, rock, np.array([0., 0.5 * dy, -0.5 * dz]))
dat.grid.add_block(bcblk)
area = 100.
con = t2connection([bcblk, blk], 1,
                   [geo.atmosphere_connection, 0.5 * dx],
                   area, 0.)
dat.grid.add_connection(con)

q0 = -0.005
gen = t2generator(block = dat.grid.blocklist[-2].name,
                  name = 'prd 1', type = 'MASS',
                  gx = q0)
dat.add_generator(gen)

dat.write(model_name + '.dat')

inc = dat.grid.incons([P0, T0])
inc.write(model_name + '.incon')

dat.run(simulator = AUTOUGH2, silent = True)

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x',
               file_format = 'gmsh22')
jsondata = dat.json(geo, mesh_filename, bdy_incons = inc,
                    mesh_coords = 'xz')
jsondata['mesh']['thickness'] = dy
jsondata['eos'] = {'name': 'se', 'initial': 'pressure'}
jsondata['thermodynamics'] = {'name': 'IAPWS'}
jsondata['initial']['primary'] = [P0, T0]
jsondata['initial']['region'] = 3
jsondata['boundaries'][0]['region'] = 3
jsondata['gravity'] = None
jsondata['output']['final'] = False
jsondata['output']['fields'] = {'fluid': ['liquidlike_fraction']}
json.dump(jsondata, open(model_name + '.json', 'w'), indent = 2, sort_keys = True)

os.chdir(orig_dir)
