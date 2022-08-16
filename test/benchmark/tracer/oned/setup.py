# 1-D single-phase and two-phase liquid tracer
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

model_name = 'oned'
dimensions = [100., 1., 1.]
nblks = [10, 1, 1]
gridsizes = [dim / nblk for dim, nblk in zip(dimensions, nblks)]
dx = [[gridsize] * nblk for gridsize, nblk in zip(gridsizes, nblks)]

geo = mulgrid().rectangular(dx[0], dx[1], dx[2], atmos_type = 2)
geo.write('g'+ model_name +'.dat')
mesh_filename = 'g' + model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x',
               file_format = 'gmsh2-binary')

phases = ['single', 'two']
initial = {'single': [30.e5, 20.], 'two': [1.e5, 0.5]}
genrate = {'single': -10. / hr, 'two': -0.1 / hr}
tstop = {'single': 100 * day, 'two': 300 * day}

dat = t2data()
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)
rock = dat.grid.rocktypelist[0]
rock.porosity = 0.1
rock.density = 2500.
rock.permeability = np.ones(3) * 1.e-13
rock.conductivity = 1.0
rock.specific_heat = 1000.

dat.relative_permeability = {'type': 1,
                             'parameters': [0., 0., 1., 1., 0., 0., 0.]}
dat.capillarity = {'type': 1,
                   'parameters': [0., 0., 1., 0., 0., 0., 0.]}

dat.lineq = {'type': 1, 'epsilon': 1.e-11}

dat.parameter.update({
     'print_level': 1,
     'relative_error': 1.e-9,
     'absolute_error': 1.
     })

dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 1 # generation table interpolation
dat.parameter['option'][24] = 0 # initial output

# inflow BC:
col = geo.columnlist[0]
lay = geo.layerlist[-1]
blkname = geo.block_name(lay.name, col.name)
blk = dat.grid.block[blkname]
bcblk = t2block('bdy 1', geo.atmosphere_volume, rock,
                np.array([0., col.centre[1], lay.centre]))
dat.grid.add_block(bcblk)
con = t2connection([bcblk, blk], 1, [geo.atmosphere_connection, 0.5 * gridsizes[0]],
                   gridsizes[1] * gridsizes[2], 0.)
dat.grid.add_connection(con)

# outflow BC:
col = geo.columnlist[-1]
lay = geo.layerlist[-1]
blkname = geo.block_name(lay.name, col.name)
gen = t2generator(name = blkname, block = blkname, type = 'MASS')
dat.add_generator(gen)

for phase in phases:

    case_model_name = '%s_%s_phase' % (model_name, phase)
    case_title = '1-D %s-phase liquid tracer problem' % phase 
    dat.title = case_title + ', steady state'
    dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
                 'num_equations': 2, 'num_secondary_parameters': 6}
    dat.generatorlist[-1].gx = genrate[phase]
    ndt = 100
    dat.parameter.update({
        'max_timesteps': ndt,
        'tstop': 1.e15,
        'print_interval': ndt,
        'default_incons': initial[phase]})
    dat.parameter['option'][16] = 5 # time step control
    inc = dat.grid.incons(initial[phase])

    jsondata = dat.json(geo, mesh_filename, incons = inc, mesh_coords = 'xz')
    jsondata['mesh']['thickness'] = dimensions[1]
    jsondata['mesh']['zones'] = {"all": {"type": "box"}}
    del jsondata['rock']['types'][0]['cells']
    jsondata['rock']['types'][0]['zones'] = 'all'
    del jsondata['output']['filename']
    jsondata['output']['initial'] = False
    jsondata['logfile'] = {'echo': False}
    json.dump(jsondata, open(case_model_name + '_ss.json', 'w'), indent = 2)
    subprocess.call(['waiwera', case_model_name + '_ss.json'])
    
    dat.simulator = 'AUTOUGH2.2EWT'
    dat.multi = {'eos': 'EWT', 'num_components': 2, 'num_phases': 2,
             'num_equations': 3, 'num_secondary_parameters': 6}
    tracer_incons = initial[phase] + [0]
    dat.parameter.update({
        'default_incons': tracer_incons})
    dat.write(case_model_name + '_ss.dat')
    inc = dat.grid.incons(tracer_incons)
    inc.write(case_model_name + '_ss.incon')
    dat.run(simulator = AUTOUGH2,
            incon_filename = case_model_name + '_ss.incon',
            silent = True)

    # transient model:
    dat.parameter.update(
        {'max_timesteps': 50,
         'tstop': tstop[phase],
         'print_interval': 1,
         'const_timestep': 10. * day,
        })
    dat.parameter['option'][16] = 0
    dat.title = case_title + ', transient'
    dat.write(case_model_name + '.dat')

    # tracer inflow:
    inflow_mass_fraction = 0.01
    inc = t2incon(case_model_name + '_ss.save')
    blk  = dat.grid.blocklist[-1]
    inc[blk.name][2] = inflow_mass_fraction
    inc.write(case_model_name + '.incon')

    dat.run(simulator = AUTOUGH2,
            incon_filename = case_model_name + '.incon',
            silent = True)

    jsondata['output']['initial'] = True
    jsondata['output']['frequency'] = dat.parameter['print_interval']
    jsondata['output']['fields'] = {'source':
                                    ['component', 'rate', 'enthalpy',
                                     'tracer_flow']}
    jsondata['time']['step']['size'] = dat.parameter['const_timestep']
    jsondata['time']['stop'] = dat.parameter['tstop']
    jsondata['time']['step']['maximum']['number'] = dat.parameter['max_timesteps']
    jsondata['time']['step']['adapt']['on'] = False

    jsondata['tracer'] = {"name": "tracer"}
    jsondata['boundaries'][0]['tracer'] = inflow_mass_fraction

    jsondata['initial'] = {'filename': case_model_name + '_ss.h5'}
    jsondata['logfile'] = {'echo': True}

    json.dump(jsondata, open(case_model_name + '.json', 'w'), indent = 2)

os.chdir(orig_dir)
