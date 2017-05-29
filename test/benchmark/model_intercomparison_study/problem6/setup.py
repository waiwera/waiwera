from t2data_json import *
from t2incons import *
from t2thermo import *
from scipy.optimize import fsolve
import json
import os

AUTOUGH2 = 'AUTOUGH2_42D'

model_dir = './run'
orig_dir = os.getcwd()
if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'problem6'

lx, ly = 5.e3, 4.e3
nx, ny = 5, 5
dx, dy = lx / float(nx), ly / float(ny)
zthick = [300.] * 4 + [600.]

geo = mulgrid().rectangular([dx] * nx, [dy] * ny, zthick, atmos_type = 1)
geo.write('g' + model_name +'.dat')

dat = t2data_export_json()
dat.title = 'Model intercomparison study problem 6'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)

rock1 = rocktype(name = 'rck 1', density = 2500., porosity = 0.2,
                permeability = np.array([100., 100., 2.]) * 1.e-15,
                specific_heat = 1000., conductivity = 1.)
dat.grid.add_rocktype(rock1)
rock2 = rocktype(name = 'rck 2', density = 2500., porosity = 0.25,
                permeability = np.array([200., 200., 50.]) * 1.e-15,
                specific_heat = 1000., conductivity = 1.)
dat.grid.add_rocktype(rock2)
for blk in dat.grid.blocklist:
    layername = geo.layer_name(blk.name)
    if layername in [' 1', ' 5']: blk.rocktype = rock1
    else: blk.rocktype = rock2
dat.grid.clean_rocktypes()

dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.1, 0., 0., 0.]}
dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}

def round4(x): return float('%.4g' % x)

ndt = 200
yr = 365. * 24. * 60 * 60.
P0, T0 = 64.e5, 280.
dt = round4(0.05 * yr)

dat.parameter.update(
    {'max_timesteps': ndt,
     'tstop': round4(6.85 * yr),
     'print_interval': 1,
     'gravity': 9.8,
     'default_incons': [P0, T0],
     'const_timestep': dt,
     'max_timestep': dt,
     })

dat.parameter['option'][1] = 1
dat.parameter['option'][11] = 2 # permeability weighting
dat.parameter['option'][12] = 1 # generation rate interpolation
dat.parameter['option'][16] = 5 # adaptive time stepping
dat.parameter['option'][24] = 2 # initial output

dat.multi = {'eos': 'EW', 'num_components': 1, 'num_phases': 2,
             'num_equations': 2, 'num_secondary_parameters': 6}

blk = geo.block_name_containing_point(np.array([500., 400., -1050.]))
gentime = [0., 2. * yr, 4. * yr, 6. * yr, 8. * yr]
genrate = [-1000., -2500., -4000., -6000.]
genrate.append(genrate[-1])
gen = t2generator(name = 'gen 1', block = blk, type = 'MASS', gx = genrate[0])
gen.time = gentime
gen.rate = genrate
gen.ltab = len(gentime)
dat.add_generator(gen)

inc = dat.grid.incons()

def layer_inc(layername, incs):
    layer = geo.layer[layername]
    for col in geo.columnlist:
        blkname = geo.block_name(layer.name, col.name)
        inc[blkname] = incs
    
T4 = 280.
P4 = sat(T4)
S4 = 0.1
d4 = cowat(T4, P4)[0]
layer_inc(' 2', [P4, S4])

T5 = 160.
def f(P): return P - P4 + 1470. * (d4 + cowat(T5, P)[0])
P5 = fsolve(f, P4)[0]
d5 = cowat(T5, P5)[0]
layer_inc(' 1', [P5, T5])

T3 = 280.
def f(P): return P - P4 - 1470. * (d4 + cowat(T3, P)[0])
P3 = fsolve(f, P4)[0]
d3 = cowat(T3, P3)[0]
layer_inc(' 3', [P3, T3])

T2 = 280.
def f(P): return P - P3 - 1470. * (d3 + cowat(T2, P)[0])
P2 = fsolve(f, P3)[0]
d2 = cowat(T2, P2)[0]
layer_inc(' 4', [P2, T2])

T1 = 280.
def f(P): return P - P2 - 1470. * (d2 + 2. * cowat(T1, P)[0])
P1 = fsolve(f, P2)[0]
d1 = cowat(T1, P1)[0]
layer_inc(' 5', [P1, T1])

Ttop = 100.
Ptop = P5 - 1470. * d5
layer_inc(' 0', [Ptop, Ttop])

# bottom BCs:
botlayer = geo.layerlist[-1]
Tbottom = 280.
Pbottom = P1 + 2940. * d1
layername = ' 6'
for col in geo.columnlist:
    blkname = geo.block_name(layername, col.name)
    blk = t2block(blkname, geo.atmosphere_volume, rock1,
                  centre = np.array(list(col.centre) + [botlayer.bottom]))
    dat.grid.add_block(blk)
    inc[blkname] = [Pbottom, Tbottom]
    botblkname = geo.block_name(botlayer.name, col.name)
    botblk = dat.grid.block[botblkname]
    con = t2connection([blk, botblk], 3,
                       [geo.atmosphere_connection, 0.5 * botlayer.thickness],
                       col.area, -1.0)
    dat.grid.add_connection(con)

# side BCs:
rightcols = [col for col in geo.columnlist if col.centre[0] > lx - dx]
for lay in geo.layerlist[1:]:
    for col in rightcols:
        edgeblkname = geo.block_name(lay.name, col.name)
        newcolname = 'Z' + col.name[1:]
        blkname = geo.block_name(lay.name, newcolname)
        if lay.name in [' 1', ' 5']: rt = rock1
        else: rt = rock2
        blk = t2block(blkname, geo.atmosphere_volume, rt,
                  centre = np.array([lx, col.centre[1], lay.centre]))
        dat.grid.add_block(blk)
        inc[blkname] = inc[edgeblkname].variable
        edgeblk = dat.grid.block[edgeblkname]
        area = dy * lay.thickness
        con = t2connection([blk, edgeblk], 1,
                           [geo.atmosphere_connection, 0.5 * dx],
                           area, 0.0)
        dat.grid.add_connection(con)

inc.write(model_name + '.incon')
dat.write(model_name + '.dat')

dat.run(simulator = AUTOUGH2,
        incon_filename = model_name + '.incon',
        silent = True)

mesh_filename = 'g' + model_name + '.exo'
geo.write_exodusii(mesh_filename)
jsondata = dat.json(geo, mesh_filename, incons = inc, bdy_incons = inc)
json.dump(jsondata, file(model_name + '.json', 'w'), indent = 2)
                  
os.chdir(orig_dir)
