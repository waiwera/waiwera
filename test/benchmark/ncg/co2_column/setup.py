# Vertical column CO2 test from O'Sullivan et al. (1985)
from __future__ import print_function
from t2data_json import *
from t2incons import *
from t2thermo import *
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'co2_column'

dx = 100.
dy = dx
dz = 20.
h = 1000.
nz = int(h / dz)
geo = mulgrid().rectangular([dx], [dy], [dz] * nz, atmos_type = 0)
geo.write('g' + model_name + '.dat')
mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')

dat = t2data_export_json()
dat.title = 'CO2 column'
dat.simulator = 'AUTOUGH2.2EWC'

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

# 'X' relative permeability curves (linear):
dat.relative_permeability = {'type': 1, 'parameters': [0.35, 0., 1., 0.7, 0.]}
dat.capillarity = {'type': 8, 'parameters': [0., 0., 0., 0., 0.]}

dat.multi = {'eos': 'EWC', 'num_components': 2, 'num_phases': 2,
             'num_equations': 3, 'num_secondary_parameters': 6}
dat.start = True

g = 9.8
ndt = 500
P0, T0 = 1.e5, 10.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e13,
     'print_interval': 1,
     'gravity': g,
     'const_timestep': 1.e5,
     'default_incons': [P0, T0, 0.]
     })
dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 2 # initial output

inc = dat.grid.incons([P0, T0, 0.])
# approximate initial conditions:
rho = 998.
for blk in dat.grid.blocklist[geo.num_atmosphere_blocks:]:
    d = -blk.centre[2]
    P = P0 + rho * g * d
    inc[blk.name] = [P, T0, 0.]

col = geo.columnlist[0]
qs = 2.0e-6 # kg/s/m2
enthalpy = 1300.e3
q0 = qs * col.area

gen = t2generator(block = dat.grid.blocklist[-1].name,
                  name = 'inj 1', type = 'MASS',
                  gx = q0, ex = enthalpy)
dat.add_generator(gen)
co2gen = t2generator(block = dat.grid.blocklist[-1].name,
                  name = 'inj 2', type = 'COM2',
                     gx = 0., ex = enthalpy)
dat.add_generator(co2gen)

CO2_mass_fractions = [0, 0.1, 1, 5] # percent

for xgp in CO2_mass_fractions:

    case_model_name = model_name + '_' + str(xgp)
    inc.write(case_model_name + '.incon')

    xg = xgp * 0.01
    gen.gx = (1. - xg) * q0
    co2gen.gx = xg * q0

    dat.write(case_model_name + '.dat')
    dat.run(simulator = AUTOUGH2, silent = True)

    jsondata = dat.json(geo, mesh_filename, mesh_coords = 'xz', incons = inc,
                        bdy_incons = inc)
    jsondata['mesh']['thickness'] = dy
    jsondata['time']['step']['solver']['linear'] = {'type': 'bcgs',
                                                    'preconditioner': {'type': 'bjacobi'}}
    jsondata['time']['step']['solver']['nonlinear']['minimum'] = {'iterations': 1} 
    json.dump(jsondata, file(case_model_name + '.json', 'w'), indent = 2, sort_keys = True)

os.chdir(orig_dir)
