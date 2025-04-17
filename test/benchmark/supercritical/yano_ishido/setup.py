# Yano-Ishido radial supercritical problem - temperature cases
from __future__ import print_function
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os
import numpy as np
from math import log
import mulgrids as mg
import IAPWS97 as I97

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'yanoish'

dr0 = 0.1
ncells = 40
alpha = 1.3
thick = 100
dr = np.logspace(start = log(dr0) / log(alpha),
                 stop = log(dr0 * alpha ** ncells) / log(alpha),
                 num = ncells, base = alpha, endpoint = False)
dy = [1.]
dz = [thick]
geo = mg.mulgrid().rectangular(dr, dy, dz)
mesh_filename = model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x', file_format = 'gmsh22')

sim = {}
sim['title'] = 'Yano-Ishido problem'
sim['mesh'] = {'filename': mesh_filename, 'radial': True}
sim['eos'] = {'name': 'se',
              'relative_permeability_modifier': 'linear',
              'conditions': 'pressure'}

# rock properties:
skin_r = 0.9
sim['mesh']['zones'] = {
    'well': {'x': [0, dr0]},
    'skin': {'x': [dr0, skin_r]},
    'reservoir': {'-': ['well', 'skin']}}
rho, hcap, cond = 2700, 1000, 2.5
sim['rock'] = {'types': [
    {'name': 'well', 'permeability': 1e-10, 'porosity': 0.99, 'zones': 'well',
     'density': rho, 'specific_heat': hcap, 'wet_conductivity': cond},
    {'name': 'skin', 'permeability': 1e-12, 'porosity': 0.05, 'zones': 'skin',
     'density': rho, 'specific_heat': hcap, 'wet_conductivity': cond},
    {'name': 'rock', 'permeability': 1e-14, 'porosity': 0.05, 'zones': 'reservoir',
     'density': rho, 'specific_heat': hcap, 'wet_conductivity': cond}],
               'relative_permeability': {
                   'type': 'linear',
                   'liquid': [0.3, 0.7],
                   'vapour': [0.3, 0.7]}
               }

# production:
sim['source'] = [{'rate': -15.7, 'cell': 0}]

# time stepping:
sim['time'] = {
    'step': {
        'size': 1,
        'adapt': {'on': True},
        'maximum': {'size': 1e5, 'number': 1000}
    },
    'stop': 1e6,
}

sim['output'] = {'initial': True, 'frequency': 1}

P0 = 30e6
for T0 in np.arange(200., 550., 50.):
    region = I97.region(T0, P0)
    sim['initial'] = {'primary': [P0, T0], 'region': region}
    if int(T0) == 350: sim['thermodynamics'] = {'extrapolate': True}
    elif 'thermodynamics' in sim: del sim['thermodynamics']
    case_model_name = '%s_%d.json' % (model_name, int(T0))
    json.dump(sim, open(case_model_name, 'w'), indent = 2, sort_keys = True)

os.chdir(orig_dir)
