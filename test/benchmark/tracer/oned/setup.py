# 1-D single-phase liquid tracer problem
from __future__ import print_function
from t2data import *
from t2incons import *
from t2thermo import cowat
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os
import subprocess

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

hr = 60. * 60.
day = 24. * hr

model_name = 'oned'
dimensions = [100., 1., 1.]
nblks = [10, 1, 1]
gridsizes = [dim / nblk for dim, nblk in zip(dimensions, nblks)]
dx = [[gridsize] * nblk for gridsize, nblk in zip(gridsizes, nblks)]

geo = mulgrid().rectangular(dx[0], dx[1], dx[2], atmos_type = 2)
geo.write('g'+ model_name +'.dat')

dat = t2data()
dat.title = '1-D single-phase liquid tracer problem'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)
rock = dat.grid.rocktypelist[0]
rock.porosity = 0.1
rock.density = 2500.
rock.permeability = np.ones(3) * 1.e-13
rock.conductivity = 1.0
rock.specific_heat = 1000.

dat.relative_permeability = {'type': 1,
                             'parameters': [0., 0., 1., 1., 0., 0., 0.]}
dat.capillarity = {'type': 1,
                   'parameters': [0., 0., 1., 0., 0., 0., 0.]}

dat.lineq = {'type': 1, 'epsilon': 1.e-11}

ndt = 100
P0, T0 = 30.e5, 20.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e15,
     'print_interval': ndt,
     'print_level': 1,
     'default_incons': [P0, T0],
     'const_timestep': 1.e6,
     'relative_error': 1.e-9,
     'absolute_error': 1.
     })

dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 1 # generation table interpolation
dat.parameter['option'][16] = 5 # time step control
dat.parameter['option'][24] = 0 # initial output

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}

# inflow BC:
col = geo.columnlist[0]
lay = geo.layerlist[-1]
blkname = geo.block_name(lay.name, col.name)
blk = dat.grid.block[blkname]
bcblk = t2block('bdy 1', geo.atmosphere_volume, rock,
                np.array([0., col.centre[1], lay.centre]))
dat.grid.add_block(bcblk)
con = t2connection([bcblk, blk], 1, [geo.atmosphere_connection, 0.5 * gridsizes[0]],
                   gridsizes[1] * gridsizes[2], 0.)
dat.grid.add_connection(con)

# outflow BC:
genrate = -10. / hr
col = geo.columnlist[-1]
lay = geo.layerlist[-1]
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = blkname, block = blkname,
                  gx = genrate, type = 'MASS')
dat.add_generator(gen)

inc = dat.grid.incons([P0, T0])

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x',
               file_format = 'gmsh2-binary')

jsondata = dat.json(geo, mesh_filename, incons = inc, mesh_coords = 'xz')
jsondata['mesh']['thickness'] = dimensions[1]
jsondata['mesh']['zones'] = {"all": {"type": "box"}}
del jsondata['rock']['types'][0]['cells']
jsondata['rock']['types'][0]['zones'] = 'all'
del jsondata['output']['filename']
jsondata['output']['initial'] = False
jsondata['logfile'] = {'echo': False}
json.dump(jsondata, file(model_name + '_ss.json', 'w'), indent = 2)

# add tracer to TOUGH2 model:
dat.simulator = 'AUTOUGH2.2EWT'
tracer_incons = [P0, T0, 0]
dat.multi = {'eos': 'EWT', 'num_components': 2, 'num_phases': 2,
             'num_equations': 3, 'num_secondary_parameters': 6}
dat.parameter.update({
    'default_incons': tracer_incons})
dat.write(model_name + '_ss.dat')

inc = dat.grid.incons(tracer_incons)
inc.write(model_name + '_ss.incon')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '_ss.incon',
        silent = True)

subprocess.call(['waiwera', model_name + '_ss.json'])

# transient model:
ndt = 5
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 100. * day,
     'print_interval': 1,
     'const_timestep': 10. * day,
     })
dat.parameter['option'][16] = 0
dat.write(model_name + '.dat')

# tracer inflow:
inflow_mass_fraction = 0.01
inc = t2incon(model_name + '_ss.save')
blk  = dat.grid.blocklist[-1]
inc[blk.name][2] = inflow_mass_fraction
inc.write(model_name + '.incon')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '.incon',
        silent = True)

jsondata['output']['initial'] = True
jsondata['output']['frequency'] = dat.parameter['print_interval']
jsondata['time']['step']['size'] = dat.parameter['const_timestep']
jsondata['time']['stop'] = dat.parameter['tstop']
jsondata['time']['step']['maximum']['number'] = ndt
jsondata['time']['step']['adapt']['on'] = False

jsondata['tracer'] = {"name": "1"}
jsondata['boundaries'][0]['tracer'] = inflow_mass_fraction

jsondata['initial'] = {'filename': model_name + '_ss.h5'}
jsondata['logfile'] = {'echo': True}

json.dump(jsondata, file(model_name + '.json', 'w'), indent = 2)

os.chdir(orig_dir)
