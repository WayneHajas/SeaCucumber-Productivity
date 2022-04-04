'''
Read data required to run the WCH productivity model.  Save it in a pickle file.

Coastlengths are hard-coded into this file through the CoastLength-dictionary.
The key-values for CoastLength also control which sites are included in the final dataset.
    i.e. if you want to exclude a site from the analysis, just remove its coastlength.

'''
import pickle
import sys
from numpy.random import seed
from datetime import datetime

sys.path.append('../../pyfunctions')
from EFAsite_moreRelBM import  Site
from EFAsite_moreRelBM import VirginSite
from constants import ConstantValues
from GetCoastLength import GetCoastLength


timestamp=datetime.now()
mdbfile='c:\StrippedBioDatabases\SeaCuke_Bio.mdb'
Survey='Zeballos'

#Coast lengths for Zeballos
CoastLength=GetCoastLength(Survey)

global ConstantValues
ConstantValues.MinYear=1901
ConstantValues.MaxYear=9999

SiteNumber=list(CoastLength.keys())
SiteNumber.sort()
Sites=[ Site(mdbfile,Survey,s,CoastLength[s])       for s in SiteNumber]
dummy=Sites[-1].GetEventTimes()
VirginSites=[ VirginSite(s)   for s in Sites]

pickle.dump((SiteNumber,Sites,VirginSites,CoastLength,timestamp),open('Zeballos.pickle','wb'))

