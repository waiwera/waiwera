# TOUGH2 radial heat pipe problem
from __future__ import print_function
from t2data import *
from t2incons import *
import json
import os
from scipy.optimize import brentq
from math import log10

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'heat_pipe'

# calculate radial mesh size expansion factors:
dr1 = 0.3
def fn(f):
    s = 0.
    fi = 1.
    for i in range(0, 100):
        s += fi
        fi *= f
    s *= dr1
    return s - 100.
f1 = brentq(fn, 1., 1.1)
dr99 = dr1 * f1 ** 99

def fn(f):
    s = 0.
    fi = f
    for i in range(0, 20):
        s += fi
        fi *= f
    return 100. + dr99 * s - 1.e4
f2 = brentq(fn, 1., 1.5)
dr100 = dr99 * f2
dr120 = dr99 * f2 ** 20

dr = np.logspace(log10(dr1), log10(dr99), num = 100)
dr = np.concatenate((dr, np.logspace(log10(dr100), log10(dr120), num = 20)))
dz = [4.5]

geo = mulgrid().rectangular(dr, [1.], dz, atmos_type = 2)
geo.write('g' + model_name + '.dat')

dat = t2data()
dat.title = '1-D radial heat pipe'
dat.simulator = 'AUTOUGH2.2EWAV'

dat.grid = t2grid().radial(dr, dz, atmos_type = 2)

rock = dat.grid.rocktypelist[0]
rock.density = 2550.
rock.porosity = 0.1
rock.permeability = np.ones(3) * 20.e-15
rock.conductivity = 2.0
rock.specific_heat = 800.

dat.relative_permeability = {'type': 7, 'parameters': [0.45, 9.6e-4, 1., 0., 0.]}
dat.capillarity = {'type': 7, 'parameters': [0.45, 1.e-3, 8.e-5, 5.e8, 1.]}

dat.multi = {'eos': 'EWAV', 'num_components': 2, 'num_phases': 2,
             'num_equations': 3, 'num_secondary_parameters': 6}
dat.start = True
dat.output_times = {'num_times_specified': 3,
                    'time': [3.15576e7, 1.2559e8, 3.15576e8]}

ndt = 250
P0, Sv0, T0 = 1.e5, 0.2, 18.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 3.15576e8,
     'print_interval': ndt,
     'gravity': None,
     'default_incons': [P0, Sv0, T0],
     'const_timestep': 1.e3,
     'relative_error': 1.e-5,
     'absolute_error': 1.,
     'derivative_increment': 1.e-7
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][16] = 4 
dat.parameter['option'][19] = 1 # alternative initial conditions
dat.parameter['option'][20] = 1 # VPL off
dat.parameter['option'][24] = 0 # initial output

# BC block at r = 10000:
r = 10.e3
col = geo.columnlist[-1]
lay = geo.layerlist[-1]
blkname = geo.block_name(lay.name, col.name)
blk = dat.grid.block[blkname]
vol = 0.
bcblk = t2block('bc %2d' % 0, vol, rock,
                np.array([r, 0., -0.5 * dz[0]]))
dat.grid.add_block(bcblk)
area = 2. * np.pi * r * dz[0]
con = t2connection([blk, bcblk], 1,
                   [0.5 * dr[-1], geo.atmosphere_connection],
                   area, 0.)
dat.grid.add_connection(con)

# Heat generator:
gen = t2generator(block = dat.grid.blocklist[0].name,
                  name = 'HTR 1',
                  type = 'HEAT',
                  gx = 3.e3)
dat.add_generator(gen)

dat.write(model_name + '.dat')
dat.run(simulator = AUTOUGH2, silent = True)
var = t2incon(model_name + '.save')[-1].variable
inc = dat.grid.incons(var)

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')
jsondata = dat.json(geo, mesh_filename, bdy_incons = inc,
                    mesh_coords = 'rz')
jsondata['initial']['primary'] = var
jsondata['initial']['region'] = 4
jsondata['mesh']['radial'] = True
jsondata['gravity'] = None
jsondata['output']['fields'] = {'fluid': ['vapour_air_mass_fraction']}
json.dump(jsondata, file(model_name + '.json', 'w'), indent = 2)

os.chdir(orig_dir)
