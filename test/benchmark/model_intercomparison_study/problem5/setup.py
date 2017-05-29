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

model_name = 'problem5'
dimensions = [300., 200., 100.]
nblks = [12, 8, 1]
gridsizes = [dim // nblk for dim, nblk in zip(dimensions, nblks)]
dx = [[gridsize] * nblk for gridsize, nblk in zip(gridsizes, nblks)]

geo = mulgrid().rectangular(dx[0], dx[1], dx[2], atmos_type = 2)
geo.write('g'+ model_name +'.dat')

dat = t2data_export_json()
dat.title = 'Model intercomparison study problem 5'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)
rock = dat.grid.rocktypelist[0]
rock.porosity = 0.35
rock.density = 2500.
rock.permeability = np.ones(3) * 25.e-15
rock.conductivity = 1.0
rock.specific_heat = 1000.

dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.1, 0., 0., 0.]}
dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}

ndt = 200
yr = 365. * 24. * 60 * 60.
P0, T0 = 36.e5, 240.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 10. * yr,
     'print_interval': 1,
     'default_incons': [P0, T0],
     'const_timestep': 0.05 * yr,
     'max_timestep': 0.05 * yr,
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

genpos = np.array([62.5, 62.5])
col = geo.column_containing_point(genpos)
genrate = -50.e-3 * dimensions[2]
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = blkname, block = blkname,
                  gx = genrate, type = 'MASS')
dat.add_generator(gen)

# BC blocks at x = 300:
rect = [np.array([dimensions[0] - gridsizes[0], 0.]), np.array([dimensions[0], dimensions[1]])]
cols = geo.columns_in_polygon(rect)
for i, col in enumerate(cols):
    blkname = geo.block_name(lay.name, col.name)
    blk = dat.grid.block[blkname]
    bcblk = t2block('bc %2d' % i, geo.atmosphere_volume, rock,
                    np.array([dimensions[0], col.centre[1], -0.5 * dimensions[2]]))
    dat.grid.add_block(bcblk)
    con = t2connection([blk, bcblk], 1, [0.5 * gridsizes[0], geo.atmosphere_connection],
                       gridsizes[1] * gridsizes[2], 0.)
    dat.grid.add_connection(con)

dat.write(model_name + 'a.dat')

# incons:
def initial_temp(pos):
    r = np.linalg.norm(pos)
    if r <= 100.: T = 240.
    elif r < 300.:
        R = (r - 100.) / 200.
        T = 240. - 160. * R * R + 80. * R**4
    else: T = 160.
    return T
inc = t2incon()
for blk in dat.grid.blocklist:
    T0 = initial_temp(blk.centre[:2])
    inc[blk.name] = [P0, T0]
inc.write(model_name + '.incon')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '.incon',
        silent = True)

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2)
jsondata = dat.json(geo, mesh_filename, incons = inc, bdy_incons = inc,
                    mesh_coords = 'xy')
jsondata['mesh']['thickness'] = dimensions[-1]
json.dump(jsondata, file(model_name + 'a.json', 'w'), indent = 2)

# problem 5b:

genpos = np.array([162.5, 137.5])
col = geo.column_containing_point(genpos)
genrate = 30.e-3 * dimensions[2]
d,u = cowat(80., P0)
h = u + P0 / d
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = blkname, block = blkname, type = 'MASS',
                  ex = h, ltab = 3,
                  time = [0., yr, dat.parameter['tstop']],
                  rate = [0., genrate, genrate])
dat.add_generator(gen)

dat.write(model_name + 'b.dat')
dat.run(simulator = 'AUTOUGH2_41Da',
        incon_filename = model_name + '.incon',
        silent = True)

jsondata = dat.json(geo, mesh_filename, incons = inc, bdy_incons = inc,
                    mesh_coords = 'xy')
jsondata['mesh']['thickness'] = dimensions[-1]
json.dump(jsondata, file(model_name + 'b.json', 'w'), indent = 2)
                  
os.chdir(orig_dir)
