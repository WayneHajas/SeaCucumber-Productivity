
import pickle
from sys import setrecursionlimit
setrecursionlimit(20000)
import sys
from numpy.random import seed
rs=20180824
seed(rs)

from pymc import  *
from pylab import *
from tables import *
import tables

project_path='//dcBCpbsNA01a/shellfish$/sf_SeaCuke$/analyses/EFA_Productivity_Model.20180929'
pyfunctions=project_path + '/pyfunctions'
case_path  =project_path + '/Jervis/NewModel_wideSdYear_WideSiteArea'
pickle_file=case_path + '/Jervis Inlet.pickle'


sys.path.append(pyfunctions)


from ExecuteModel import ExecuteModel
from EFAsite_moreRelBM import Site
from EFAProdModel_wideSdYear_WideSiteArea import EFAProdModel

SiteNumber,Sites,VirginSites,CoastLength,timestamp=pickle.load(open(pickle_file,"rb"))

Model=EFAProdModel(SiteNumber,Sites,VirginSites,CoastLength,oldB=1)
print('\nModel is Made\n')

Niter=210000
burn=10000
thin=10


name=case_path + '/seed.'+str(rs)

db='hdf5'
verbose=2
maxSimplex=0
maxPowell=0
maxBFGS=0

ExecuteModel(Model,Niter=Niter,burn=burn,thin=thin,name=name,db=db,verbose=2,maxSimplex=maxSimplex,maxPowell=maxPowell,maxBFGS=maxBFGS)
