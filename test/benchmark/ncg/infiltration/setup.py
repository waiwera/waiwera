# TOUGH2 1-D infiltration problem
from __future__ import print_function
from t2data_json import *
from t2incons import *
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

AUTOUGH2 = 'AUTOUGH2_42D'
model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'infiltration'

dx, nx = 0.005, 40
geo = mulgrid().rectangular([dx] * nx, [1.], [1.], atmos_type = 2)
mulgrid_format_specification['node'][1][1:] = ['10.3f'] * 2
geo.write('g' + model_name + '.dat')

dat = t2data_export_json()
dat.title = '1-D infiltration'
dat.simulator = 'AUTOUGH2.2EWAV'

dat.grid = t2grid().fromgeo(geo)

rock = dat.grid.rocktypelist[0]
rock.density = 2385.
rock.porosity = 0.45
rock.permeability = np.ones(3) * 1.2e-14
rock.conductivity = 1.045
rock.specific_heat = 1030.

dat.relative_permeability = {'type': 1, 'parameters': [0.333, -0.1, 1., 0., 0.]}
dat.capillarity = {'type': 1, 'parameters': [9.7902e3, 0.333, 1., 0., 0.]}

dat.multi = {'eos': 'EWAV', 'num_components': 2, 'num_phases': 2,
             'num_equations': 3, 'num_secondary_parameters': 6}
dat.start = True
dat.output_times = {'num_times_specified': 3, 'time': [864., 5184., 9504.]}

ndt = 100
P0, Sv0, T0 = 1.e5, 0.56, 20.
Pa0 = .97663438450445e5

dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e4,
     'print_interval': ndt,
     'gravity': None,
     'default_incons': [P0, Sv0, T0],
     'const_timestep': 100.,
     'max_timestep': 500.,
     'relative_error': 1.e-5,
     'absolute_error': 1.,
     'derivative_increment': 1.e-7
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][16] = 4 
dat.parameter['option'][19] = 1 # alternative initial conditions
dat.parameter['option'][20] = 1 # VPL off
dat.parameter['option'][24] = 2 # initial output

# BC block at r = 0:
col = geo.columnlist[0]
lay = geo.layerlist[-1]
blkname = geo.block_name(lay.name, col.name)
blk = dat.grid.block[blkname]
vol = 0.
bc_name = 'bc  0'
bcblk = t2block(bc_name, vol, rock, np.array([0., 0.5, -0.5]))
dat.grid.add_block(bcblk)
area = 1.
con = t2connection([bcblk, blk], 1,
                   [geo.atmosphere_connection, 0.5 * dx],
                   area, 0.)
dat.grid.add_connection(con)

dat.write(model_name + '.dat')

inc = dat.grid.incons([P0, Sv0, T0])
inc[bc_name] = [P0, T0, 0.]
inc.write(model_name + '.incon')

dat.run(simulator = AUTOUGH2, silent = True)

# Waiwera input:
mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')
jsondata = dat.json(geo, mesh_filename, bdy_incons = inc,
                    mesh_coords = 'xz')
jsondata['initial']['primary'] = [P0, Sv0, Pa0]
jsondata['initial']['region'] = 4
jsondata['boundaries'][0]['region'] = 1
jsondata['gravity'] = None
jsondata['output']['final'] = False
json.dump(jsondata, file(model_name + '.json', 'w'), indent = 2, sort_keys = True)

os.chdir(orig_dir)
