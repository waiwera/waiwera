# 1-D tracer doublet problem
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
yr = 365.25 * day

model_name = 'doublet'
dimensions = [200., 10., 10.]
nblks = [100, 1, 1]
gridsizes = [dim / nblk for dim, nblk in zip(dimensions, nblks)]
dx = [[gridsize] * nblk for gridsize, nblk in zip(gridsizes, nblks)]

geo = mulgrid().rectangular(dx[0], dx[1], dx[2], atmos_type = 2)
geo.write('g'+ model_name +'.dat')

dat = t2data()
dat.title = '1-D tracer doublet problem'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)
rock = dat.grid.rocktypelist[0]
rock.porosity = 0.1
rock.density = 2600.
rock.permeability = np.array([2.e-13, 2.e-15, 2.e-15])
rock.conductivity = 1.5
rock.specific_heat = 900.

dat.relative_permeability = {'type': 1,
                             'parameters': [0.5, 0., 1., 0.5, 0., 0., 0.]}
dat.capillarity = {'type': 1,
                   'parameters': [0., 0., 1., 0., 0., 0., 0.]}

dat.lineq = {'type': 1, 'epsilon': 1.e-11}

ndt = 500
P0, T0 = 50.e5, 10.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e15,
     'print_interval': ndt,
     'print_level': 1,
     'default_incons': [P0, T0],
     'const_timestep': 1.e4
     })

dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 2 # generation table interpolation
dat.parameter['option'][16] = 5 # time step control
dat.parameter['option'][24] = 0 # initial output

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}

lay = geo.layerlist[-1]

# injection:
q = 0.5
h = 4.21e5
col = geo.columnlist[0]
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = 'inj 1', block = blkname,
                  type = 'MASS', gx = q, ex = h)
dat.add_generator(gen)

# production:
PI = 4.e-11
Pc = 50.e5
qmax = 5.0
col = geo.columnlist[-1]
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = 'prd 1', block = blkname,
                  type = 'DELT', gx = PI, ex = Pc, hg = qmax)
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
json.dump(jsondata, open(model_name + '_ss.json', 'w'), indent = 2)

# add tracer to TOUGH2 model:
dat.simulator = 'AUTOUGH2.2EWTD'
tracer_incons = [P0, T0, 0]
dat.multi = {'eos': 'EWTD', 'num_components': 2, 'num_phases': 2,
             'num_equations': 3, 'num_secondary_parameters': 8}
dat.parameter.update({
    'default_incons': tracer_incons})
D = 1.e-4
dat.diffusion = [[-D, -D], [-D, -D]]
dat.write(model_name + '_ss.dat')

inc = dat.grid.incons(tracer_incons)
inc.write(model_name + '_ss.incon')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '_ss.incon',
        silent = True)
t2incon(model_name + '_ss.save').write(model_name + '.incon')

subprocess.call(['waiwera', model_name + '_ss.json'])

# transient model:
ndt = 300
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 0.5 * yr,
     'print_interval': 10,
     'const_timestep': 0.01 * day,
     'max_timestep': 1 * day,
     })

inflow_mass_fraction = 1.e-3
end_time = 0.15 * day
qt = q * inflow_mass_fraction
col = geo.columnlist[0]
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = 'trc 1', block = blkname,
                  type = 'COM2')
gen.ltab = 4
gen.itab = '1'
gen.time = [0., end_time, end_time + 1.e-3, 1. * yr]
gen.rate = [qt, qt, 0., 0.]
gen.enthalpy = [h] * 4
dat.add_generator(gen)

dat.write(model_name + '.dat')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '.incon',
        silent = True)

jsondata['output']['frequency'] = dat.parameter['print_interval']
jsondata['output']['fields'] = {'source':
                                 ['component', 'rate', 'enthalpy',
                                  'tracer1_flow']}
jsondata['time']['step']['size'] = dat.parameter['const_timestep']
jsondata['time']['stop'] = dat.parameter['tstop']
jsondata['time']['step']['maximum']['number'] = ndt
jsondata['time']['step']['maximum']['size'] = dat.parameter['max_timestep']

jsondata['tracer'] = {'name': 'tracer1', "diffusion": D}
jsondata['source'][0]['tracer'] = [[0, qt], [end_time, 0]]
jsondata['source'][0]['interpolation'] = 'step'

jsondata['initial'] = {'filename': model_name + '_ss.h5'}
jsondata['logfile'] = {'echo': True}

json.dump(jsondata, open(model_name + '.json', 'w'), indent = 2)

os.chdir(orig_dir)
