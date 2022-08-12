# Waiwera demo problem

import json
import layermesh.mesh as lm

# set up mesh:
dx = [1000] * 10
dy = dx
dz = [10] * 5 + [20] * 5 + [50] * 5 + [100] * 5 + [200] * 5
mesh = lm.mesh(rectangular = (dx, dy, dz))
mesh.translate((0,0,200))
surface_data = [(0,0,0), (10e3,0,0), (10e3,10e3,0), (0,10e3,0,0),
                (5e3,5e3,200)]
mesh.fit_surface(surface_data)
mesh.write('demo_mesh.h5')
mesh_filename = 'demo_mesh.exo'
mesh.export(mesh_filename)

sim = {}
sim['title'] = 'Waiwera demo problem'
sim['mesh'] = {'filename': mesh_filename}

# initial and boundary conditions:
P0, T0 = 1.e5, 20.
sim['initial'] = {'primary': [P0, T0]}
sim['boundaries'] = [{'primary': [P0, T0], 'region': 1,
                    'faces': {'cells': [c.index for c in mesh.surface_cells],
                              'normal': [0, 0, 1]}}]

# rock properties:
sim['mesh']['zones'] = {'cap': {'z': [0, 200]}, 'reservoir': {'-': 'cap'}}
sim['rock'] = {'types': [
    {'name': 'rock', 'permeability': 20.e-15, 'zones': 'reservoir'},
    {'name': 'caprock', 'permeability': 0.5e-15, 'zones': 'cap'}
]}

# upflow:
upflow_cols = mesh.find([(4000,4000), (6000,6000)])
upflow_cells = [col.cell[-1].index for col in upflow_cols]
sim['source'] = [{'rate': 100, 'enthalpy': 1300e3,
                  'cells': upflow_cells}]

# time stepping:
sim['time'] = {
    'step': {
        'size': 1e6,
        'adapt': {'on': True},
        'maximum': {'number': 100}
    },
    'stop': 1e16,
}

sim['output'] = {'initial': False, 'frequency': 0}

json.dump(sim, open('demo.json', 'w'), indent = 2, sort_keys = True)
