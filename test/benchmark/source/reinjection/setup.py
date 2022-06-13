# 1-D column problem with reinjection

from t2data import *
from t2incons import *
from t2thermo import sat, cowat
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os

AUTOUGH2 = 'AUTOUGH2'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'reinjection'

d = 100.
dx = [d]
dy = [d]
dz = list(range(40, 140, 10)) + [150]

geo = mulgrid().rectangular(dx, dy, dz, atmos_type = 0)
geo.write('g'+ model_name +'.dat')
mesh_filename = 'g' + model_name + '.exo'
geo.write_mesh(mesh_filename)

title = '1-D column problem with reinjection'
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
dat.multi = {'num_components':1, 'num_equations':2,
             'num_phases':2, 'num_secondary_parameters':6,
             'eos': 'EW'}

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
ndt = 100
day = 24. * 60. * 60.
dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': 365 * day,
     'const_timestep': 0.5 * day,
     'max_timestep': 7. * day,
     'print_interval': 1
    })

# production wells:
production_wellnames = ['prd 1', 'prd 2', 'prd 3']
wells = {'prd 1': (-180, 0.55e6), 'prd 2': (-250, 0.6e6), 'prd 3': (-350, -1)}
PI = 5e-12
Pb = 10e5
qsmax = 0.
col = geo.columnlist[0]
for wellname in production_wellnames:
    z, Psep = wells[wellname]
    lay = geo.layer_containing_elevation(z)
    blkname = geo.block_name(lay.name, col.name)
    prd = t2generator(name = wellname, block = blkname,
                      gx = PI, ex = Pb, fg = Psep, hg = qsmax,
                      type = 'DELG')
    dat.add_generator(prd)

jsondata = dat.json(geo, mesh_filename, incons = inc)
jsondata['output']['fields'] = {'source': ['natural_cell_index',
                                           'rate', 'enthalpy',
                                           'steam_rate'],
                                'network_group': ['rate', 'enthalpy',
                                                  'steam_rate']}

# setup source network manually (until PyTOUGH can convert it):
for g in jsondata['source'][1: 1 + len(production_wellnames)]:
    Psep = wells[g['name']][1]
    if Psep < 0: Psep = [1.45e6, 0.55e6]
    g['separator'] = {'pressure': Psep}

# reinjection:
z = -20
hs = 85.e3
lay = geo.layer_containing_elevation(z)
blkname = geo.block_name(lay.name, col.name)
cell_index = geo.block_name_index[blkname] - geo.num_atmosphere_blocks

qf = 1.5
inj = t2generator(name = 'inj 1', block = blkname,
                  gx = qf, ex = hs, hg = -1., type = 'FINJ')
dat.add_generator(inj)
jsondata['source'].append({'name': 'inj 1', 'cell': cell_index})

pf = 0.1
inj = t2generator(name = 'inj 2', block = blkname,
                  ex = hs, hg = -pf, type = 'PINJ')
dat.add_generator(inj)
jsondata['source'].append({'name': 'inj 2', 'cell': cell_index})

z = -65
lay = geo.layer_containing_elevation(z)
blkname = geo.block_name(lay.name, col.name)
cell_index = geo.block_name_index[blkname] - geo.num_atmosphere_blocks
hl = 440.e3
qmax = 2.0
Pref = 9.e5
injectivity = 3.e-6
inj = t2generator(name = 'inj 3', block = blkname,
                  gx = qmax, ex = hl, hg = Pref, fg = injectivity, type = 'IMAK')
dat.add_generator(inj)
jsondata['source'].append({'name': 'inj 3', 'cell': cell_index,
                           'direction': 'injection',
                           'enthalpy': hl,
                           'limiter': {'total': qmax},
                           'injectivity': {'pressure': Pref,
                                           'coefficient': injectivity}})
overflow_pf = 0.05
inj = t2generator(name = 'inj 4', block = blkname,
                  ex = hl, hg = overflow_pf, fg = 1., type = 'RINJ')
dat.add_generator(inj)
jsondata['source'].append({'name': 'inj 4', 'cell': cell_index})

jsondata['network'] = {
    'group': [
        {'name': 'reinjection group',
         'in': production_wellnames}],
    'reinject': [
        {'name': 'reinjector', 'in': 'reinjection group', 'overflow': 'overflow reinjector',
         'water': [
             {'out': 'inj 3'}
         ],
         'steam': [
             {'out': 'inj 1', 'rate': qf, 'enthalpy': hs},
             {'out': 'inj 2', 'proportion': pf, 'enthalpy': hs}
         ]},
        {'name': 'overflow reinjector',
         'water': [{'out': 'inj 4', 'proportion': overflow_pf, 'enthalpy': hl}]
        }
    ]}

case_name = 'run'
dat.title = title
jsondata['title'] = dat.title
jsondata['output']['filename'] = model_name + '.h5'
dat.write(model_name + '.dat')
inc.write(model_name + '.incon')
dat.run(simulator = AUTOUGH2, silent = True)
json.dump(jsondata, open(model_name + '.json', 'w'),
          indent = 2, sort_keys = True)
