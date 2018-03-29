from t2data_json import *
from t2incons import *
from builtins import zip
from t2thermo import cowat
from scipy.optimize import fsolve
import json
import os

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'problem4'
depth = 2000.
num_cells = 20
width = 100.
dz = depth / num_cells
g = 9.8

geo = mulgrid().rectangular([width], [width], [dz] * num_cells, atmos_type = 0)
geo.write('g' + model_name + '.dat')

dat = t2data_export_json()
dat.title = 'Model intercomparison study problem 4'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)
rock = dat.grid.rocktypelist[0]
rock.name = 'upper'
rock.porosity = 0.15
rock.permeability = np.ones(3) * 5.e-15
rock.density = 2500.
rock.conductivity = 1.0
rock.specific_heat = 1000.

rock = rocktype()
rock.name = 'lower'
rock.porosity = 0.25
rock.permeability = np.ones(3) * 100.e-15
rock.density = 2500.
rock.conductivity = 1.0
rock.specific_heat = 1000.
dat.grid.add_rocktype(rock)
for blk in dat.grid.blocklist[1:]:
    if blk.centre[2] < -1000.:
        blk.rocktype = dat.grid.rocktype['lower']

dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.05, 0., 0., 0.]}
dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}

year = 365. * 24. * 60. * 60.
P0, T0 = 1.013e5, 10.
dat.parameter.update(
    {'max_timesteps': 150,
     'tstop': 40. * year,
     'max_timestep': 0.5 * year,
     'print_interval': 1,
     'gravity': g,
     'default_incons': [P0, T0],
     'const_timestep': 1.e5,
     'relative_error': 1.e-6,
     'absolute_error': 1.
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 1 # generation table interpolation
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 0 # initial output

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}

lay = geo.layerlist[-1]
col = geo.columnlist[0]
genrate = -100. / 1.e6 * col.area
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = blkname, block = blkname,
                  gx = genrate, type = 'MASS')
dat.add_generator(gen)

dat.write(model_name + '.dat')

# incons:
inc = t2incon()
def temp(d):
    if d <= 1000.: return 10. + 280. * d / 1000.
    else: return 270. + 20. * d / 1000.
depth, T, P = [0.], [T0], [P0]
rho = cowat(T0, P0)[0]
depth.append(0.5 * dz)
T.append(temp(depth[-1]))         
def f(p): return p - (P0 + g * 0.5 * (rho + cowat(T[-1], p)[0]) * depth[-1])
P.append(fsolve(f, P0)[0])
for i in range(num_cells):
    rho = cowat(T[-1], P[-1])[0]
    depth.append(depth[-1] + dz)
    T.append(temp(depth[-1]))
    def f(p): return p - (P[-1] + g * 0.5 * (rho + cowat(T[-1], p)[0]) * dz)
    P.append(fsolve(f, P[-1])[0])
    
for blk, Tblk, Pblk in zip(dat.grid.blocklist, T, P):
    inc[blk.name] = [Pblk, Tblk]
inc.write(model_name + '.incon')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '.incon',
        silent = True)

mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')
jsondata = dat.json(geo, mesh_filename, incons = inc, bdy_incons = inc,
                    mesh_coords = 'xz')
jsondata['mesh']['thickness'] = width
jsondata['output']['fields'] = {'fluid': ['liquid_saturation']}
json.dump(jsondata, file(model_name + '.json', 'w'), indent = 2)

os.chdir(orig_dir)
