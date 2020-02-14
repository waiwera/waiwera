# 1-D MINC column problem

from t2data import *
from t2incons import *
from t2thermo import sat, cowat
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'minc_column'

d = 100.
dx = [d]
dy = [d]
dz = list(range(40, 140, 10)) + [150]

geo = mulgrid().rectangular(dx, dy, dz, atmos_type = 0)
geo.write('g'+ model_name +'.dat')
mesh_filename = 'g' + model_name + '.exo'
geo.write_mesh(mesh_filename)

title = '1-D MINC column problem'
dat = t2data()
dat.title = title + ': steady state'
dat.simulator = 'AUTOUGH2.2EW'
dat.start = True
dat.grid = t2grid().fromgeo(geo)

dat.grid.rocktype['dfalt'].permeability = np.ones(3) * 1.e-13

ndt = 50
P0, T0 = 1.e5, 20.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e15,
     'print_interval': ndt,
     'const_timestep': 1.e6,
     'gravity': 9.8,
     'relative_error': 1.e-5,
     'timestep_reduction': 4.0,
     'default_incons': [P0, T0],
     'max_iterations': 8
     })
dat.parameter['option'][1] = 1
dat.parameter['option'][16] = 6
dat.parameter['option'][24] = 2
dat.start = False

dat.relative_permeability = {'type': 1, 'parameters': [0.,0.,1.,1.]}
dat.capillarity = {'type': 1, 'parameters': [0.,0.,0.,0.]}

gx = 10. # kg/s
t = 240. # deg C
p = sat(t)
d,u = cowat(t,p)
h = u + p / d
gen = t2generator(name = 'gen 1', block = dat.grid.blocklist[-1].name,
                  gx = gx, ex = h, type = 'MASS')
dat.add_generator(gen)

# steady state run:
dat.write(model_name + '_ss.dat')
inc = dat.grid.incons([P0, T0])
inc.write(model_name + '_ss.incon')
dat.run(simulator = AUTOUGH2, silent = True)
inc = t2incon(model_name + '_ss.save')

# transient runs:
ndt = 40
day = 24. * 60. * 60.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 90 * day,
     'const_timestep': 0.5 * day,
     'max_timestep': 3. * day,
     'print_interval': 1
    })

z= -350.
gx = -25.
lay = geo.layer_containing_elevation(z)
col = geo.columnlist[0]
blkname = geo.block_name(lay.name, col.name)
prd = t2generator(name = 'prd 1', block = blkname,
                  gx = gx, type = 'MASS')
dat.add_generator(prd)

case_name = 'single'
dat.title = title + ': single porosity'
dat.write(model_name + '_' + case_name + '.dat')
inc.write(model_name + '_' + case_name + '.incon')
dat.run(simulator = AUTOUGH2, silent = True)
jsondata = dat.json(geo, mesh_filename, incons = inc)
json.dump(jsondata, file(model_name + '_' + case_name + '.json', 'w'),
          indent = 2, sort_keys = True)

vol = [0.1, 0.3, 0.6]
nf = 3
fracture_rockname, matrix_rockname = 'fract', 'matrx'
def minc_rockname(rockname, level):
    return matrix_rockname if level > 0 else fracture_rockname

spacing = 5
minc_zone_bounds = [-600, -100]
minc_layers = [lay for lay in geo.layerlist
               if minc_zone_bounds[0] < lay.centre < minc_zone_bounds[1]]
minc_blks = [geo.block_name(lay.name, col.name) for lay in minc_layers]
minc_level, minc_inc = dat.grid.minc(vol, spacing, nf,
                                     minc_rockname = minc_rockname, blocks = minc_blks,
                                     incon = inc)
dat.grid.clean_rocktypes()
case_name = 'minc'
dat.title = title + ': MINC'

rock = dat.grid.rocktype['dfalt']
fracture = dat.grid.rocktype[fracture_rockname]
fracture.porosity = 0.5
matrix = dat.grid.rocktype[matrix_rockname]
matrix.porosity = (rock.porosity - fracture.porosity * vol[0]) / (1. - vol[0])
matrix.permeability = np.ones(3) * 1.e-18

dat.write(model_name + '_' + case_name + '.dat')
minc_inc.write(model_name + '_' + case_name + '.incon')
dat.run(simulator = AUTOUGH2, silent = True)

jsondata = dat.json(geo, mesh_filename, incons = inc)
jsondata['mesh']['zones'] = {'minc': {'z': minc_zone_bounds}}
jsondata['mesh']['minc'] = {
    'geometry': {
        'fracture': {'volume': vol[0], 'planes': nf, 'spacing': spacing},
        'matrix': {'volume': vol[1:]}},
    'rock': {'fracture': {'type': fracture_rockname},
             'matrix': {'type': matrix_rockname},
             'zones': 'minc'}}

json.dump(jsondata, file(model_name + '_' + case_name + '.json', 'w'),
          indent = 2, sort_keys = True)
