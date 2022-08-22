# One-cell tracer decay
from __future__ import print_function
import mulgrids
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')
import os
import subprocess

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

hr = 60. * 60.
day = 24. * hr

model_name = 'decay'
h = 10.

geo = mulgrids.mulgrid().rectangular([h], [h], [h], atmos_type = 2)
mesh_filename = model_name + '.msh'
geo.write_mesh(mesh_filename, dimension = 2, slice = 'x',
               file_format = 'gmsh2-binary')

P = 10.e5
T = 60.
region = 1
X = 0.001
k0 = 1.e-6
Ea = 2.e3

inp = {}
inp['title'] = 'one-cell tracer decay'
inp['mesh'] = {'filename': mesh_filename, 'thickness': h}
inp['eos'] = 'we'
inp['initial'] = {'primary': [P, T], 'region': region, 'tracer': X}
inp['rock'] = {
    'relative_permeability': {'type': 'linear'},
    'types': [{'name': 'rock', 'porosity': 0.1, 'cells': [0]}]
}
inp['output'] = {'initial': True, 'frequency': 1, 'final': True}

inp['time'] = {
    'step': {
        'method': 'bdf2',
        'size': 1 * day,
        'adapt': {'on': False},
        'maximum': {'number': 20},
        },
    'stop': 20 * day
    }

inp['tracer'] = [
    {'name': 'no_decay'},
    {'name': 'constant', 'decay': k0},
    {'name': 'temperature', 'decay': k0, 'activation': Ea}
    ]

json.dump(inp, open(model_name + '.json', 'w'), indent = 2)

os.chdir(orig_dir)
