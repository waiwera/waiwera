# Vertical column salt CO2 test
from __future__ import print_function
from t2data import *
from t2incons import *
from t2thermo import *
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

AUTOUGH2 = 'AUTOUGH2'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'salt_co2_column'

dx = 100.
dy = dx
dz = [30.] * 10 + [35.] * 20
geo = mulgrid().rectangular([dx], [dy], dz, atmos_type = 0)
geo.write('g' + model_name + '.dat')
mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x', file_format = 'gmsh22')

dat = t2data()
dat.title = 'Salt CO2 column'
dat.simulator = 'AUTOUGH2.2EWSG'

dat.grid = t2grid().fromgeo(geo)

rock = dat.grid.rocktypelist[0]
rock.permeability = np.ones(3) * 20.e-15
rock.conductivity = 2.0

caprock = rocktype(name = 'capr',
                   permeability = np.ones(3) * 0.5e-15,
                   conductivity = 2.0)
dat.grid.add_rocktype(caprock)
for blk in dat.grid.blocklist[1:]:
    if blk.centre[2] > -300.: blk.rocktype = caprock

dat.relative_permeability = {'type': 1, 'parameters': [0.35, 0., 1., 0.7, 0.]}
dat.capillarity = {'type': 8, 'parameters': [0., 0., 0., 0., 0.]}

dat.multi = {'eos': 'EWSG', 'num_components': 3, 'num_phases': 3,
             'num_equations': 4, 'num_secondary_parameters': 6}
dat.selection['integer'] = [0] * 9 + [0] + [0] + [0] * 5 + [2]
dat.start = True

g = 9.8
ndt = 200
P0, T0 = 1.e5, 20.
Xs0, Xg0 = 0., 0.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e15,
     'print_interval': 1,
     'gravity': g,
     'const_timestep': 1.e5,
     'default_incons': [P0, Xs0, Xg0, T0]
     })
dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 2 # initial output

inc = dat.grid.incons([P0, Xs0, Xg0, T0])
# approximate initial conditions:
rho = 998.
for blk in dat.grid.blocklist[geo.num_atmosphere_blocks:]:
    d = -blk.centre[2]
    P = P0 + rho * g * d
    inc[blk.name] = [P, Xs0, Xg0, T0]

col = geo.columnlist[0]
qs = 2.0e-6 # kg/s/m2
enthalpy = 1800.e3
q0 = qs * col.area
xs = 0.1
xg = 0.05

gen = t2generator(block = dat.grid.blocklist[-1].name,
                  name = 'inj 1', type = 'MASS',
                  gx = (1 - xs - xg) * q0, ex = enthalpy)
dat.add_generator(gen)
salt_gen = t2generator(block = dat.grid.blocklist[-1].name,
                       name = 'inj 2', type = 'COM2',
                       gx = xs * q0, ex = enthalpy)
dat.add_generator(salt_gen)
co2gen = t2generator(block = dat.grid.blocklist[-1].name,
                  name = 'inj 3', type = 'COM3',
                     gx = xg * q0, ex = enthalpy)
dat.add_generator(co2gen)

case_model_name = model_name
inc.write(case_model_name + '.incon')

dat.write(case_model_name + '.dat')
dat.run(simulator = AUTOUGH2, silent = True)

dat.simulator = 'AUTOUGH2.2EW'
dat.multi['eos'] = 'EW'
jsondata = dat.json(geo, mesh_filename, mesh_coords = 'xz', incons = inc)
jsondata['mesh']['thickness'] = dy
jsondata['time']['step']['solver']['linear'] = {'type': 'bcgs',
                                                'preconditioner': {'type': 'asm'}}
jsondata['time']['step']['solver']['nonlinear']['minimum'] = {'iterations': 1}
jsondata['output'] = {'fields': {'fluid': ['liquid_saturation', 'liquid_CO2_mass_fraction']}}

jsondata['eos'] = 'wsce'
jsondata['initial']['primary'] = [P0, T0, Xs0, Xg0]
jsondata['initial']['region'] = 1
jsondata['boundaries'][0]['primary'] = [P0, T0, Xs0, Xg0]
jsondata['boundaries'][0]['region'] = 1

json.dump(jsondata, open(case_model_name + '.json', 'w'), indent = 2, sort_keys = True)

os.chdir(orig_dir)
