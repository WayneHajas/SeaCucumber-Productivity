'''2016-01-06
       Force 
         *Site Effects to sum to zero
         *Year Effects to sum to zero
         *Transect Effects to sum to zero for each site
2016-01-12
    Use pymc.Deterministic to explicitly incorporate the relative-biomasses as nodes.
2016-01-13
    Get all the nodes defined in __init__    
    I would like to break up __init__ but some nodes do not seem to get sampled properly.
2016-02-15
   Make the prior distributions more like the old WinBugs version
   F:\Archive\s-drive\analyses\CukeExpHarvRecov\Zeballos.beta\ZeballosSpatial.odc for example
   Corrected serious error in calculation of virgin biomass.  Now use the zeroth year-effect
2016-02-15
    Start from the new implementatino of the WCh model and create teh Gram-schaeffer model
2016-02-29
	Apply a prior specifying that virgin grand-mean must be approximate to the nonvirgin grand-mean.
	Apply a prior to sdYear to try and make it small
 2016-03-04
     Implemented virgin transect-counts as random variables
2016-03-29
    VPredVal:  make conversion between population and biomass 
2017-01-17
    Change the name of the model to EFAProdModel and start a new implementation
    
2017-02-06
    Change constructor-variables so be SiteNumber,Sites,VirginSites instead of mdbfile,Survey,CoastLength
    
    Previous versions read from a .mdb file - and that required a 32-bit version of python.
    This version assumes recieves data in the form of instances of classes.  No need to read .mdb files.  Can use the 64-bit version of python.
2017-03-13
    Do not create year-effects for years where there was harvest but no survey    
2017-03-24
    Apply prior-distributions to relative-biomass values.  
    dummy~normal(log(RBM),0.754)   
    
2017-04-07
    Replaced in error in cacluating b from binit and xmax
2017-04-11
    Error Correction.  Relative biomass included in calculation of predicted number of animals on transect
2017-09-13
    Just Refactoring code
2017-10-18
    Create a variable and method, SiteAreasFromData that will sample site-area from data.
    This value can be compared against estimates that are influenced by harvest/productivity/other surveys
    
2017-11-07
    simplify EFAProdModel.SiteAreas
    Rather than a prior-distribution based on a transformed beta-variable, just
    use a uniform distribution with range of +/-two standard errors
2017-11-09
    Correct implementation of prior for last element in an array of effects.
    Use LastEffect distribution
    Need a new, nonsensical, distribution for the dummy data.  
    Probability gets larger away from mean.
2018-01-05
    Correct VirginBiomass
    Year Effects do not affect virgin biomass.
	
2018-01-08
	Adjust priors on productivity parameters and on effect standard deviations to be the same as used in WinBUGS version
2018-06-18
    Change priors on MeanWeight.  Use truncated normal.  Bounds correspond to 95% confidence bounds from field-data.
2019-06-03
	Special version with wide distribution of sdyear
2020-06-23
    noSiteAreaBoundary
        Prior distribution for mean-transect-length is simple normal.  mu and sigma based on length of virgin-transects
2020-06-29
         Prior distribution of mean transect-length is based on +/-3 standard errors
		 Use  modified version of EFsite to save post-harvest values of RelBM
'''

from pymc import  *
from pylab import *
from tables import *
from numpy import log
import tables
import warnings
from datetime import datetime as dt

from ADO import adoBaseClass as daoBaseClass
from constants import ConstantValues
from EFAsite_moreRelBM import  Site,VirginSite
from Harvest import  Harvest
from LastEffect import LastEffect
from MakeName import MakeName
from MeanWeight import MeanWeight
from NextRelBmass import CurRelBmass
from NodeNameToDecimalYear import MatchNodeRefTime
from ToDecimalYear import ToDecimalYear,FromDecimalYear

from EFAProdModel import EFAProdModel as oldEFAProdModel




class EFAProdModel(oldEFAProdModel):
    def __init__(self,SiteNumber,Sites,VirginSites,CoastLength,oldB=1):
        '''EFAProdModel(SiteNumber,Sites,VirginSites,CoastLength,oldB=1)
        * Sites and VirginSites are lists of instances of the Site and VirginSite classes
        * CoastLength is a dictionary giving lengths in metres for each site. e.g.{0:10000,2:10002,4:10004,8:1008,16:10016}            
        * oldB is the original relative-abundance.   Value is probably one
        
         '''
        self.SiteNumber=SiteNumber
        self.nsite=len(self.SiteNumber)
        self.Sites=Sites
        self.VirginSites=VirginSites
        self.CoastLength=CoastLength
        self.oldB=oldB


        self.ProductivityParameters()  #Generate nodes for productivity function       
        self.GrandMeans()              #Generate nodes for Grand Mean and virgin Grand Mean of biomass density
        self.EffectStandardDeviations()#Generate nodes for standard deviations of effects
        self.SiteEffects()             #                   site effects
        self.TransectEffects()         #                   transect effects
        self.YearEffects()             #                   Year Effects        
        self.SiteAreas()               #                   Site Areas
        self.VirginBiomass()           #                   Virgin biomass
        self.CalcRelAbund()            #                   Relative Abundance
                              
        # Mean Weight nodes will be created as they are needed
        self.MeanWeight =[]
             
        self.VirginTransects()#Nodes to represent expected and observed number of animals for initial survey (untrimmed transects)
        self.SurveyTransects()#Nodes to represent expected and observed number of animals for trimmed transects
        self.SiteAreasFromData()#Nodes to represent site-areas as sampled directly from statistics
        

        print('Finished Constructing Model')
    

    def EffectStandardDeviations(self):
        #nodes representing variance
        self.sdYear=  Lognormal('sdYear',    -1.38,0.52)
        self.sdSite=    Uniform('sdSite',    .01,2.3)
        self.sdTransect=Uniform('sdTransect', .01,2.3)
        self.tauYear    =Lambda('tauYear',     lambda x=self.sdYear:    1./x/x)
        self.tauSite    =Lambda('tauSite',     lambda x=self.sdSite:    1./x/x)
        self.tauTransect=Lambda('tauTransect', lambda x=self.sdTransect:1./x/x)
    

    def SiteAreas(self): 
        #Site areas
        nsite=len(self.SiteNumber)
        w=[ self.CoastLength[s]  for s in self.SiteNumber]
        l=[ t.AvgTranLen  for t in self.VirginSites]
        se=[ t.sterrTranLen  for t in self.VirginSites]
		#Limit range of mean transect-length to +/- three standard errofs
        self.AvgTran=[Uniform( MakeName('AvgTran',self.SiteNumber[i]), (l[i]-3*se[i]), (l[i]+3*se[i])  ) for i in range(nsite)]
		
        self.sArea=[Lambda(MakeName('sArea',self.SiteNumber[i]),lambda \
            a=self.AvgTran[i],w=w[i]: a*w) for i in range(nsite)]
        
      
    def SiteAreasFromData(self): 
        #Site areas
        nsite=len(self.SiteNumber)
        w=[ self.CoastLength[s]  for s in self.SiteNumber]
        l=[ t.AvgTranLen  for t in self.VirginSites]
        se=[ t.sterrTranLen  for t in self.VirginSites]
        tau=[1/s/s  for s in se]
                
        self.AvgTranLenFromData=[ Normal(MakeName('AvgTranLenFromData',self.Sites[i].SiteNumber),l[i],tau[i]  )   for i in range(nsite)]
        self.sAreaFromData=[Lambda(MakeName('sAreaFromData',self.SiteNumber[i]),lambda \
            a=self.AvgTranLenFromData[i],w=w[i]: a*w) for i in range(nsite)]
        
        
         
        
        

if __name__ == "__main__":
  import pickle
  pfile="../MCMC/Jervis/Jervis.pickle"
  SiteNumber,Sites,VirginSites,CoastLength=pickle.load(open(pfile,"rb"))  
  
  ConstantValues.MinYear=1999  #Earliest year to be considered in calculations
  
  test=EFAProdModel(SiteNumber,Sites,VirginSites,CoastLength,oldB=1)
  #print(test.GetMeanWeight(16,[2005,2,14]))  

