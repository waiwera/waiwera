# Supercritical pure water vertical column problem
from t2data import *
from t2incons import *
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

AUTOUGH2 = 'AUTOUGH2_2S'
model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

h = 3500.e3 # upflow enthalpy, kJ/kg
q = 1.e-7 # upflow mass flow rate kg/m2/s
k = 0.1 # permeability mD
P0, T0 = 1.e5, 20.

model_name = 'colv'
dz = np.logspace(1, 2.6, 50)
dx = 10.
dy = dx
geo = mulgrid().rectangular([dx], [dy], dz, atmos_type = 0, convention = 1)
geo.write('g' + model_name + '.dat')

dat = t2data()
dat.simulator = 'AUTOUGH2.2EW'
dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}
dat.start = True

dat.grid = t2grid().fromgeo(geo)
rock = dat.grid.rocktypelist[0]
rock.permeability = np.ones(3) * k * 1.e-15

dat.relative_permeability = {'type': 1, 'parameters': [0.3, 0., 1., 0.7, 0.]}
dat.capillarity = {'type': 8, 'parameters': [0., 0., 0., 0., 0.]}

inc = dat.grid.incons([P0, T0])

g = 9.8
ndt = 1000
fdinc = 1.e-9
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 1.e15,
     'print_interval': ndt,
     'gravity': g,
     'const_timestep': 1.e6,
     'derivative_increment': fdinc,
     'print_block': dat.grid.blocklist[-1].name
     })
dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][16] = 5
dat.parameter['option'][24] = 2 # initial output
dat.parameter['default_incons'] = [P0, T0]

blkname = dat.grid.blocklist[-1].name
gen = t2generator(name = blkname, block = blkname, type = 'MASS',
                  gx = geo.columnlist[0].area * q, ex = h)
dat.add_generator(gen)

dat.write(model_name + '.dat')
inc.write(model_name + '.incon')

mesh_filename = model_name + '.msh'
geo.write_mesh(mesh_filename, file_format = 'gmsh22')
jsondata = dat.json(geo, mesh_filename, eos = 1, incons = inc)

jsondata['thermodynamics'] = {'name': 'IAPWS'}
jsondata['eos'] = {'name': 'se'}
jsondata['initial'] = {'filename': model_name + '_cold.h5'}
jsondata['time']['step']['solver']['nonlinear']['jacobian'] = {'differencing': {'increment':
                                                                                fdinc}}

json.dump(jsondata, open(model_name + '.json', 'w'), indent = 2, sort_keys = True)

# EOS3 version
dat.simulator = 'AUTOUGH2.2EWA'
dat.multi['eos'] = 'EWA'
dat.multi['num_components'] = 2
dat.multi['num_equations'] = 3
Sv0 = 1.e-6
dat.parameter['default_incons'] = [P0, Sv0, T0]
inc = t2incon(model_name + '_cold.incon')
dat.write(model_name + '.dat')
inc.write(model_name + '.incon')
dat.run(simulator = AUTOUGH2, silent = True)

os.chdir(orig_dir)
