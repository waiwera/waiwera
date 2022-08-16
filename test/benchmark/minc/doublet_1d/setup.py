# 1-D MINC doublet problem

from t2data import *
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

yr = 365.25 * 24. * 3600.

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'minc_1d'

l = 500.
ndx = 10
delx = l / ndx
dx = [delx] * ndx
dy = [delx]
dz = [delx]
geo = mulgrid().rectangular(dx, dy, dz)
geo.write('g' + model_name + '.dat')
mesh_filename = 'g' + model_name + '.exo'
geo.write_mesh(mesh_filename)

dat = t2data()
dat.title = '1-D MINC problem'
dat.simulator = 'AUTOUGH2.2EW'
dat.start = True
dat.grid = t2grid().fromgeo(geo)

rock = rocktype('rock ', density = 2650., porosity = 0.1, 
                 permeability = np.ones(3) * 6.e-15,
                 conductivity = 2.1, specific_heat = 1000.)
dat.grid.add_rocktype(rock)
for blk in dat.grid.blocklist: blk.rocktype = rock
dat.grid.clean_rocktypes()

dat.parameter.update({
        'print_level': 1, 'print_interval': 99,
        'max_timesteps': 99, 'tstart': 0.0, 'tstop': 50. * yr,
        'const_timestep': 1.e5, 'max_timestep': 2 * yr,
        'print_block': 'AA 11', 'default_incons': [85.e5, 0.01],
        'relative_error': 1.e-5, 'derivative_increment': 1.e-8
        })
dat.parameter['option'][1] = 1
dat.parameter['option'][16] = 4

dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.05, 0., 0., 0.]}
dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}

dat.output_times = {'num_times_specified': 2, 'num_times': 2,
                    'time': [5. * yr, 25. * yr]}

layer = geo.layerlist[-1]
rate = 0.1
col = geo.columnlist[0]
blkname = geo.block_name(layer.name, col.name)
inj = t2generator(block = blkname, name = 'inj 1', type = 'MASS',
                  gx = rate, ex = 0.5e6)
dat.add_generator(inj)

col = geo.columnlist[-1]
blkname = geo.block_name(layer.name, col.name)
proj = t2generator(block = blkname, name = 'pro 1', type = 'MASS',
                   gx = -rate)
dat.add_generator(proj)

single_model_name = model_name + '_single'
dat.write(single_model_name + '.dat')
dat.run(simulator = AUTOUGH2, silent = True)

jsondata = dat.json(geo, mesh_filename)
jsondata['initial']['region'] = 4
json.dump(jsondata, open(single_model_name + '.json', 'w'),
          indent = 2, sort_keys = True)

# MINC runs:

vol = [0.1, 0.9]
nf = 3
fracture_rockname, matrix_rockname = 'fract', 'matrx'
def minc_rockname(rockname, level):
    return matrix_rockname if level > 0 else fracture_rockname

for spacing in [50, 100, 200]:

    dat = t2data(single_model_name + '.dat')
    dat.grid.minc(vol, spacing, nf, minc_rockname = minc_rockname)
    dat.grid.clean_rocktypes()
    dat.title += ': fracture spacing = %d' % spacing

    fracture = dat.grid.rocktype[fracture_rockname]
    fracture.porosity = 0.5
    matrix = dat.grid.rocktype[matrix_rockname]
    matrix.porosity = (rock.porosity - fracture.porosity * vol[0]) / (1. - vol[0])
    matrix.permeability = np.ones(3) * 1.e-18

    case_model_name = model_name + '_' + str(spacing)
    dat.write(case_model_name + '.dat')
    dat.run(simulator = AUTOUGH2, silent = True)

    jsondata = dat.json(geo, mesh_filename)
    jsondata['initial']['region'] = 4
    jsondata['mesh']['zones'] = {'all': {'-': None}}
    jsondata['mesh']['minc'] = {
        'geometry': {
            'fracture': {'volume': vol[0], 'planes': nf, 'spacing': spacing},
            'matrix': {'volume': vol[1:]}},
        'rock': {'fracture': {'type': fracture_rockname},
                 'matrix': {'type': matrix_rockname},
                 'zones': 'all'}}

    json.dump(jsondata, open(case_model_name + '.json', 'w'),
              indent = 2, sort_keys = True)
    
