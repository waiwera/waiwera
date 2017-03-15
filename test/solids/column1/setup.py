from t2data_json import *
from t2incons import *
from t2thermo import cowat
import json
import os

model_dir = './run'
orig_dir = os.getcwd()

if not os.path.isdir(model_dir): os.makedirs(model_dir)
os.chdir(model_dir)

model_name = 'column1'
dimensions = [10.,10.,30.]
nblks = [1, 1, 3]
gridsizes = [dim/nblk for dim,nblk in zip(dimensions,nblks)]
dx = [[gridsize]*nblk for gridsize,nblk in zip(gridsizes,nblks)]

geo = mulgrid().rectangular(dx[0],dx[1],dx[2],atmos_type=2)
geo.write('g'+model_name+'.dat')

dat = t2data_export_json()
dat.title = 'first solid mechanics simulation'
dat.simulator = 'AUTOUGH2.2'

dat.grid = t2grid().fromgeo(geo)
rock = dat.grid.rocktypelist[0]
rock.porosity = 0.35
rock.density = 2500.
rock.permeability = np.ones(3)*25.e-15
rock.conductivity = 1.0
rock.specific_heat = 1000.

P0 = 101325.
T0 = 20.

os.chdir(orig_dir)
