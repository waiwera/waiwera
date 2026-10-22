# Supercritical air/water 2D slice problem
from t2data import *
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

AUTOUGH2 = 'AUTOUGH2_2S'
model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

qup = 11 # upflow heat flux W/m2
k = 1.e-15 # permeability m2
kup = np.array([100., 100., 10.]) * 1.e-15 # shallow permeability m2
background = 40e-3 # background heat flux W/m2
height = 200
qx = 6.e3

model_name = 'slice'

dx = np.logspace(2, 3.25, 10)
w = np.sum(dx)
dy = 1e3
dz = np.logspace(1, 3.15, 26)
geo = mulgrid().rectangular(dx, [dy], dz, atmos_type = 1, convention = 2)

def topo(x):
    return -x * height / w

for col in geo.columnlist:
    col.surface = topo(col.centre[0])
    geo.set_column_num_layers(col)
geo.snap_columns_to_nearest_layers()
geo.setup_block_name_index()
geo.setup_block_connection_name_index()
geo.write('g' + model_name + '.dat')

dat = t2data()
dat.simulator = 'AUTOUGH2.2EWA'
dat.multi = {'eos': 'EWA', 'num_components': 2, 'num_phases': 2,
             'num_equations': 3, 'num_secondary_parameters': 6}
dat.start = True

dat.grid = t2grid().fromgeo(geo)

rock = dat.grid.rocktypelist[0]
rock.permeability = np.ones(3) * k
upper = rocktype(name = 'upper', permeability = kup)
dat.grid.add_rocktype(upper)
for blk in dat.grid.blocklist:
    if blk.centre[2] >= -300: blk.rocktype = upper

dat.relative_permeability = {'type': 1, 'parameters': [0., 0., 1., 1., 0.]}
dat.capillarity = {'type': 8, 'parameters': [0., 0., 0., 0., 0.]}

g = 9.8
ndt = 1200
P0, Sv0, T0 = 1.01325e5, 0.99, 20.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e16,
     'print_interval': ndt,
     'gravity': g,
     'const_timestep': 1.e5,
     'print_block': dat.grid.blocklist[-1].name,
     'default_incons': [P0, Sv0 + 10, T0]
     })
dat.parameter['option'][1] = 1
dat.parameter['option'][5] = 3 # region transition output
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 0 # initial output

dat.parameter['default_incons'] = [P0, 0, T0]
# inc = dat.grid.incons(dat.parameter['default_incons'])
# for blk in dat.grid.blocklist[:geo.num_atmosphere_blocks]:
#     inc[blk.name] = [P0, Sv0 + 10, T0]
inc = t2incon(model_name + '_cold.incon')

layer = geo.layerlist[-1]
for col in geo.columnlist:
    xi = col.centre[0] / qx
    if xi <= 1:
        qval = (1 - xi) * qup + xi * background
    else:
        qval = background
    blkname = geo.block_name(layer.name, col.name)
    gen = t2generator(name = blkname, block = blkname, type = 'HEAT',
                      gx = col.area * qval)
    dat.add_generator(gen)

yr = 365.25 * 24 * 3600
rainfall = 1.0
rain_density = 998.
infiltration = 0.1
hrain = 83.93e3
qrain = infiltration * rainfall * rain_density / yr
for col in geo.columnlist:
    layer_index = geo.num_layers - col.num_layers
    layer = geo.layerlist[layer_index]
    blkname = geo.block_name(layer.name, col.name)
    gen = t2generator(name = blkname, block = blkname, type = 'MASS',
                      gx = col.area * qrain,
                      ex = hrain)
    dat.add_generator(gen)

dat.write(model_name + '.dat')
inc.write(model_name + '.incon')

dat.run(simulator = AUTOUGH2, silent = True)

mesh_filename = 'g' + model_name + '.exo'
geo.write_mesh(mesh_filename)
jsondata = dat.json(geo, mesh_filename, incons = inc)
jsondata['thermodynamics'] = {'name': 'iapws'}
jsondata['eos'] = {'name': 'sae', 'primary': {'scale': {'partial_pressure': 1e5}}}
jsondata['initial'] = {'filename': model_name + '_cold.h5'}
json.dump(jsondata, open(model_name + '.json', 'w'), indent = 2, sort_keys = True)

os.chdir(orig_dir)
