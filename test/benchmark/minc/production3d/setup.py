# 3-D MINC production problem with well on deliverability, and
# time-dependent productivity index

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

model_name = 'minc_3d'

yr = 365. * 24. * 60 * 60.

dx = np.array([500, 750, 1000, 1250, 1500])
dy = np.array([400, 600, 800, 1000, 1200])
dz = np.array([350, 500, 650, 700, 800])

geo = mulgrid().rectangular(dx, dy, dz, atmos_type = 1, convention = 0)
geo.write('g' + model_name +'.dat')
mesh_filename = 'g' + model_name + '.exo'
geo.write_mesh(mesh_filename)

dat = t2data()
dat.title = '3-D MINC production model'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)

lay1 = rocktype(name = 'lay 1', porosity = 0.2, density = 2500.,
                 permeability = np.array([2., 2., 1.]) * 1.e-15,
                 conductivity = 1., specific_heat = 1000.)
dat.grid.add_rocktype(lay1)
lay2 = rocktype(name = 'lay 2', porosity = 0.25, density = 2500.,
                 permeability = np.array([10., 10., 2.]) * 1.e-15,
                 conductivity = 1., specific_heat = 1000.)
dat.grid.add_rocktype(lay2)
cap1 = rocktype(name = 'cap 1', porosity = 0.2, density = 2500.,
                 permeability = np.array([0.2, 0.2, 1.]) * 1.e-15,
                 conductivity = 1., specific_heat = 1000.)
dat.grid.add_rocktype(cap1)
cap2 = rocktype(name = 'cap 2', porosity = 0.2, density = 2500.,
                 permeability = np.array([0.2, 0.2, 2.]) * 1.e-15,
                 conductivity = 1., specific_heat = 1000.)
dat.grid.add_rocktype(cap2)

for blk in dat.grid.blocklist:
    x, y, z = blk.centre[0], blk.centre[1], blk.centre[2]
    if z < -1500: blk.rocktype = lay1
    elif z < -350.: blk.rocktype = lay2
    elif x < 1250. and y < 1000: blk.rocktype = cap2
    else: blk.rocktype = cap1
dat.grid.clean_rocktypes()

dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.1, 0., 0., 0.]}
dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}

ndt = 1000
P0, T0 = 1.e5, 15.

dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e16,
     'print_interval': ndt,
     'gravity': 9.8,
     'default_incons': [P0, T0],
     'const_timestep': 1.e6,
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 1 # generation rate interpolation
dat.parameter['option'][16] = 5 # adaptive time stepping
dat.parameter['option'][24] = 2 # initial output

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}

dat.lineq = {'type': 2, 'epsilon': 1.e-11}

# bottom inflows:

def inflow_rate(pos):
    x,y = pos[0], pos[1]
    if x < 1250 and y < 1000: q = 2.e-5
    elif x < 2250 and y < 1800: q = 5.e-6
    else: q = 4.e-8
    return q

h = 1250.e3
layer = geo.layerlist[-1]
for col in geo.columnlist:
    blkname = geo.block_name(layer.name, col.name)
    genname = geo.block_name('99', col.name)
    q = inflow_rate(col.centre) * col.area
    gen = t2generator(name = genname, block = blkname,
                      type = 'MASS', gx = q, ex = h)
    dat.add_generator(gen)

dat.write(model_name + '_ss.dat')

inc = dat.grid.incons((P0, T0))
inc.write(model_name + '_ss.incon')
dat.run(simulator = AUTOUGH2, silent = True)
inc = t2incon(model_name + '_ss.save')

# production run:

blk = geo.block_name_containing_point(np.array([10, 10, -1000.]))
gen = t2generator(name = 'gen 1', block = blk, type = 'DELG',
                  gx = 1.e-10, ex = 0.1e5, hg = 150.)
time = [0, 2 * yr, 3 * yr, 4 * yr]
PI = [1.e-11, 2.e-11, 3.e-11, 3.e-11]
gen.time = time
gen.rate = PI
gen.ltab = len(time)
dat.add_generator(gen)

def round4(x): return float('%.4g' % x)
dt = round4(0.05 * yr)
ndt = 100
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': round4(4 * yr),
     'print_interval': 1,
     'const_timestep': dt,
     'max_timestep': dt,
     })

# MINC:

vol = [0.1, 0.3, 0.6]
nf = 3
fracture_rockname, matrix_rockname = 'fract', 'matrx'
def minc_rockname(rockname, level):
    return matrix_rockname if level > 0 else fracture_rockname

spacing = 10
xmax, ymax = 2250, 1800
zrange = [-1500, -350]
minc_blks = [blk.name for blk in dat.grid.blocklist[geo.num_atmosphere_blocks:]
             if blk.centre[0] < xmax and blk.centre[1] < ymax
             and zrange[0] < blk.centre[2] < zrange[1]]
minc_level, minc_inc = dat.grid.minc(vol, spacing, nf,
                                     minc_rockname = minc_rockname, blocks = minc_blks,
                                     incon = inc)
rock = dat.grid.rocktype['lay 2']
fracture = dat.grid.rocktype[fracture_rockname]
fracture.porosity = 0.7
matrix = dat.grid.rocktype[matrix_rockname]
matrix.porosity = (rock.porosity - fracture.porosity * vol[0]) / (1. - vol[0])
matrix.permeability = np.ones(3) * 1.e-18
dat.grid.clean_rocktypes()

dat.write(model_name + '.dat')
minc_inc.write(model_name + '.incon')
dat.run(simulator = AUTOUGH2, silent = True)

jsondata = dat.json(geo, mesh_filename, incons = inc)
jsondata['mesh']['zones'] = {'minc': {'x': [0, xmax], 'y': [0, ymax], 'z': zrange}}
jsondata['mesh']['minc'] = {
    'geometry': {
        'fracture': {'volume': vol[0], 'planes': nf, 'spacing': spacing},
        'matrix': {'volume': vol[1:]}},
    'rock': {'fracture': {'type': fracture_rockname},
             'matrix': {'type': matrix_rockname},
             'zones': 'minc'}}

json.dump(jsondata, file(model_name + '.json', 'w'),
          indent = 2, sort_keys = True)
