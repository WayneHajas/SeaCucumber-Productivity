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
from EFAsite import  Site,VirginSite
from Harvest import  Harvest
from LastEffect import LastEffect
from MakeName import MakeName
from MeanWeight import MeanWeight
from NextRelBmass import CurRelBmass
from NodeNameToDecimalYear import MatchNodeRefTime
from ToDecimalYear import ToDecimalYear,FromDecimalYear




class EFAProdModel():
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
    

    def GetYearNode(self,YearValue):
        '''EFAProdModel.GetYearNode(YearValue)
        Get the year-node corresponding to a specific date
        '''
        
        ydate=FromDecimalYear(YearValue) #force to datetime-value
        if isinstance(ydate,(list,ndarray)):
            result=[self.GetYearNode(y)  for y in ydate]
            return(result)
        y=ydate.year #Calendar year
        try:
            SearchName=MakeName('year',y) #Name of node we are looking for
            result=[t for t in self.YearNode if t.__name__==SearchName][0] 
            return(result)
        except:
            #Appropriate year-node was not found
            print('\n EFAProdModel 238')
            print('No Node found for year ',YearValue)
            return(None)
    
    def GetSiteNode(self,SiteValue):
        '''EFAProdModel.GetSiteNode(SiteValue)
        Get the site-node corresponding to a site-number
        '''
        #if SiteValue is an instance of Site-class, just extract the site-number from there.
        if isinstance(SiteValue,Site):
            result=self.GetSiteNode(SiteValue.Site)
            return(result)
        try:
            SearchName=MakeName('site',SiteValue)#Name of node we are looking for
            result=[t for t in self.SiteNode if t.__name__==SearchName][0] 
            return(result)
        except:
            #Appropriate year-node was not found
            print('\n EFAProdModel 267')
            print('No Node found for site ',SiteValue)
            return(None)            
                 
    def GetVBmass(self,SiteNumber):
        '''EFAProdModel.GetVBmass(SiteNumber)
        Get the node corresponding to virgin-biomass for a site-number
        '''
        SearchName=MakeName('VBmass',SiteNumber)#Name of node we are looking for
        try:
            result=[t for t in self.VBmass if t.__name__==SearchName][0] 
            return(result)
        except:
            #Appropriate node was not found
            print('\n EFAProdModel 272')
            print('No Node found for site ',SiteValue)
            return(None)  
         
        
        
    def GetSite(self,SiteNumber):
        '''EFAProdModel.GetSite(SiteNumber)
        Get the instance of Site-class corresponding to a site-number
        '''
        match=[s for s in self.Sites if s.SiteNumber==SiteNumber ]
        if (match==[]) or (match==None):
            print('\n TotalEFAWCH 224')
            print('Did not find Site Number - ',SiteNumber)
            print('Available SiteNumbers', self.SiteNumber)
            for s in self.Sites:print(s.SiteNumber)
            return(None)
        if(len(match)>1):
            print('\n TotalEFAWCH 228')
            print('Multiple matches to Site Number - ',SiteNumber)
            return(None)
        #Everything worked OK
        return(match[0])    
    
    def GetVSite(self,SiteNumber):
        '''EFAProdModel.GetVSite(SiteNumber)
        Get the instance of VirginSite-class corresponding to a site-number
        '''
        match=[s for s in self.VirginSites if s.SiteNumber==SiteNumber ]
        if (match==[]) or (match==None):
            print('\n EFAProdModel 302')
            print('Did not find Virgin Site Number - ',SiteNumber)
            return(None)
        if(len(match)>1):
            print('\n EFAProdModel 306')
            print('Multiple matches to VirginSite Number - ',SiteNumber)
            return(None)
        #Everything worked OK
        return(match[0])    
    
    def GetHarvest(self,SiteNumber):
        '''EFAProdModel.GetHarvest(SiteNumber)
        Get the instance of harvest-class corresponding to a site-number
        '''
        CurSite=self.GetSite(SiteNumber)
        result=CurSite.SiteHarv
        return(result)
    
    def GetRelAbund(self,SiteNumber,SurveyDate):
        '''EFAProdModel.GetRelAbund(SiteNumber,SurveyDate)
        Get the relative-abundance node corresponding to a site-number and date
        '''
        
        #Construct the name the node should have.
        dttime=FromDecimalYear(SurveyDate)
        SearchName=MakeName('Rel_Abund',[SiteNumber,dttime])
        try:
            #Search in the list of relative-abundance nodes
            result=[t for t in self.RelAbund if t.__name__==SearchName][0]
            return(result)
        except:
            print('\n EFAProdModel 335')
            print('No Relative Abundance found for ',SearchName)
            print('  Site',SiteNumber)
            print('  SurveyDate',SurveyDate)
            #for t in self.RelAbund: print(t.__name__)
            return(None)
    
    def GetMeanWeight(self,SiteNumber,SurveyDate):
        '''
        EFAProdModel.GetMeanWeight(SiteNumber,SurveyDate)
        Get a mean-weight node corresponding to a site-number and a survey-date
        
        If such a node already exists, that is what is returned.  
        If the node does not exist, it is created and added to the list of available mean-weight-nodes.
        '''
        CurSite=self.GetSite(SiteNumber)
        
        #Best weight data corresponding to Survey Date.  Might not exactly match the day.
        BestData=CurSite.GetBestAvgWeight(SurveyDate)
        
        #See if there is a node corresponding to this day of best weight-data
        SearchName=MakeName('MeanWeight',[SiteNumber ,BestData['year'],BestData['month'],BestData['day'] ])
        CurNode=[ t   for t in self.MeanWeight if t.__name__==SearchName]   
        
        #Node already exists for the date and site-number
        if CurNode:
            return(CurNode[0])
            
        #Need to make a new node
        mu=BestData['AvgWgt']
        tau=BestData['tau']
        a=BestData['95CB'][0]
        b=BestData['95CB'][1]
        NewNode=TruncatedNormal(SearchName,mu,tau,a=a,b=b)
        self.MeanWeight+=[NewNode]
        return(NewNode)
            
    def GetTranEff(self,s,TransectNumber):
        '''EFAProdModel.GetTranEff(s,TransectNumber)
        Get the transect-node corresponding to a site-number and transect-number
        '''
        SearchName=MakeName('TranEffect',[s,TransectNumber]) 
        result=[TE for TE in self.TranEffect if TE.__name__==SearchName  ][0]
        return(result)
   
    def GetPredVal(self,s,TranNum,SurveyDate):
        '''EFAProdModel.GetPredVal(s,TranNum,SurveyDate)
        Get the node for predicted relative biomass corresponding to a site-number and transect-number and survey-date
        '''
        dttime=FromDecimalYear(SurveyDate)
        SearchName=MakeName('PredVal',[s,TranNum,dttime])
        result=[TE for TE in self.PredVal if TE.__name__==SearchName  ][0]
        return(result)

    def GetYearOfSurvey(self):
        '''EFAProdModel.GetYearOfSurvey()
        Get all year-values where a survey takes place
        '''
        YearSurvey=[]
        for s in self.Sites:
            YearSurvey+=s.GetSurveyYears()
        result=list(set(YearSurvey))
        result.sort()
        return(result)
        
    def ProductivityParameters(self):
        #parameters of WCH productivity function        
        self.xmax=Uniform('xmax',0.001, 0.999)        
        self.binit=Uniform('binit',0.001, 0.999)
        self.b=Lambda('b',lambda binit=self.binit,xmax=self.xmax:(xmax<=0.5)*binit + (xmax>0.5)*binit*(1-xmax)/xmax)
        self.a=Lambda('a',lambda b=self.b,xmax=self.xmax:b*xmax/(1-xmax))
        self.fmax=Uniform('fmax',.01,0.25)    
    
    def GrandMeans(self):
        #Grand mean.  animals per metre-squared.  Applies to truncated transects        
        self.lnGrandMean=Uniform('lnGrandMean',-6,5)
        self.VGrandMean=Lognormal('VGrandMean',self.lnGrandMean,10) #Applied to virgin transects (no truncation)
        
    def EffectStandardDeviations(self):
        #nodes representing variance
        self.sdYear=    Uniform('sdYear',    .01,0.5)
        self.sdSite=    Uniform('sdSite',    .01,2.3)
        self.sdTransect=Uniform('sdTransect', .01,2.3)
        self.tauYear    =Lambda('tauYear',     lambda x=self.sdYear:    1./x/x)
        self.tauSite    =Lambda('tauSite',     lambda x=self.sdSite:    1./x/x)
        self.tauTransect=Lambda('tauTransect', lambda x=self.sdTransect:1./x/x)
        
    def SiteEffects(self):
        #Generate Site Effects
        self.SiteNode=[Normal(MakeName('site',s),0,self.tauSite) for s in self.SiteNumber[:-1]    ]
        s=self.SiteNumber[-1]
        self.SiteNode+=[Lambda(MakeName('site',s),lambda x=self.SiteNode:-sum(x) ) ]#Site Effects sum to zero
        self.DummySite=LastEffect('DummySite',value=0,mu=self.SiteNode[-1],sigma=self.sdSite,n=self.nsite)#zero-valued dummy-value used to impose prior distribution on final node in list
         
       
    def YearEffects(self):
        #Generate Year Effects
        self.Year=self.GetYearOfSurvey()
        self.YearNode=[Normal(MakeName('year',s),0,self.tauYear) for s in self.Year[:-1]    ]
        s=self.Year[-1]
        self.YearNode+=[Lambda(MakeName('year',s),lambda x=self.YearNode:-sum(x))  ]#Year Effects sum to zero
        self.DummyYear=LastEffect('DummyYear',value=0,mu=self.YearNode[-1],sigma=self.sdYear,n=len(self.YearNode))#zero-valued dummy-value used to impose prior distribution on final node in list
           

    def SiteAreas(self): 
        #Site areas
        nsite=len(self.SiteNumber)
        w=[ self.CoastLength[s]  for s in self.SiteNumber]
        l=[ t.AvgTranLen  for t in self.VirginSites]
        se=[ t.sterrTranLen  for t in self.VirginSites]
        cb95=[[w[i]*(l[i]-2*se[i]), w[i]*(l[i]+2*se[i])  ] for i in range(nsite)]
        
        #Site-areas have uniform distribution  based upon +/- two standard errors
        self.sArea=[Uniform(MakeName('sArea',self.SiteNumber[i]) ,cb95[i][0],cb95[i][1]   )  for i in range(nsite)]
        
      
      
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
        
        
 
    def VirginBiomass(self):
        #Virgin biomass
        YE=[self.YearNode[0] for s in self.SiteNumber]#Same year-node used for all virgin sites
        SE=[self.GetSiteNode(s)  for s in self.SiteNumber]#Get the site-nodes
        self.VBmass=[Lambda(MakeName('VBmass',self.SiteNumber[i]), lambda \
                        GrandMean=self.VGrandMean,\
                        SiteEff=SE[i],\
                        area=self.sArea[i]:\
                        area*GrandMean*exp(SiteEff)) for i in range(self.nsite)]
 
                        
    def CalcRelAbund(self):
        #All the necessary relative abundance values.SpatiaL.  Fraction of virgin-value
        self.RelAbund=[]
        self.DummyRBM1=[]
        self.DummyRBM2=[]
        for s in self.SiteNumber:
            print('EFAProdModel 131',s)
            CurSite=self.GetSite(s)
            VBiomass=self.GetVBmass(s)#Node representing virgin biomass
            EventTimes=CurSite.GetEventTimes()
            ntime=len(EventTimes)
                     
 
            #The remainining time-values        
            for i in range(ntime):
                NewName=MakeName('Rel_Abund',[s,FromDecimalYear(EventTimes[i])]) #MakeName does not deal with real-values
                print('EFAProdModel 141 ',NewName)
                nt=EventTimes[i] #New Time in the interval
                self.RelAbund+=[\
                              pymc.Deterministic(eval=CurRelBmass,\
                              name=NewName,\
                              parents={   'NewTime':nt,
                                          'a':self.a,
                                          'b':self.b,
                                          'xmax':self.xmax,
                                          'fmax':self.fmax,
                                          'VBiomass':VBiomass,
                                          'HarvestHist':CurSite.SiteHarv.hdata
                              },
                              doc='Relative Biomass'+NewName,
                              trace=True,
                              verbose=1,
                              dtype=float,
                              plot=False,
                              cache_depth=0
                              )]  

    def TransectEffects(self):
        self.TranEffect=[] #Transect Effects   
        self.DummyTran=[]  #Used to impose prior distributions on last transect of site 
        for s in self.SiteNumber:
            SE=self.GetSiteNode(s) #Site Effect
            curSite=self.GetSite(s)
            
            #Generate transect-effects for the site 
            CurTran=[ Normal(MakeName('TranEffect',[s,t]),0,self.tauTransect)  for t in curSite.TranNum[:-1]] #All but the last transect in the site
            t=curSite.TranNum[-1]  #Last transect          
            CurTran+=[Lambda(MakeName('TranEffect',[s,t]),lambda x=CurTran:-sum(x))  ]#Transect Effects sum to zero for each site          
            self.DummyTran+=[LastEffect(MakeName('DummyTran',s),value=0,mu=CurTran[-1],sigma=self.sdTransect,n=len(CurTran))]#zero-valued dummy-value used to impose prior distributionon final node in list
            self.TranEffect+=  CurTran       

    def VirginTransects(self):
        #Virgin transects, no trimming  
     
        self.VPredVal=[]   #Predicted number of animals on a  virgin transect (no truncation)
        self.VObsVal=[]    #Observed number of animals on a virgin transect    
        for s in self.SiteNumber:
            SE=self.GetSiteNode(s) #Site Effect
            curSite=self.GetSite(s)
            YE=self.GetYearNode(curSite.DayOfSurvey[0])
            vMW=self.GetMeanWeight(s,curSite.DayOfSurvey[0]) #Mean Weights at first survey
            CurPredVirgin=[Lambda(MakeName('PredVirg',[s,vt.TranNum]),lambda area=vt.GetArea(),GM=self.VGrandMean,vMW=vMW,YE=YE,SE=SE,TE=self.GetTranEff(s,vt.TranNum):\
                    area*GM/vMW*exp(YE+SE+TE),trace=False)   for vt in curSite.VirgTran]  
            self.VPredVal+=CurPredVirgin
            CurVObsVal=[Poisson(MakeName('VObsVal',vt.__name__[8:]),vt,observed=True,value=curSite.VirgTran[i].GetNumCuke()) for i,vt in  enumerate(CurPredVirgin)]
            self.VObsVal+=CurVObsVal    

    def SurveyTransects(self):
        self.PredVal=[]    #Predicted number of animals on a transect
        self.ObsVal=[]     #Observed number of animals on a transect       
        for s in self.SiteNumber:
            SE=self.GetSiteNode(s) #Site Effect
            curSite=self.GetSite(s)
            for curTran in curSite.SurvTran:
                TranNum=curTran.TransectNumber
                TE=self.GetTranEff(s,TranNum)#Current Transect Effect
                
                useArea=curTran.GetArea()
                #Each day that the transect was surveyd                
                for DayTran in curTran.AllQuad:
                    y,m,d=DayTran.year,DayTran.month,DayTran.day
                    CurDate=dt(y,m,d) #Standard date-format
                    YE=self.GetYearNode(CurDate)
                    CurNumAnimals=DayTran.GetNumCuke()
                    
                    # suffix to node-names                    
                    subscript=  MakeName('_',[s,TranNum,y,m,d])
                    
                    #Predicted number of animals to be in the truncated-transect  
                    RBM=self.GetRelAbund(s,CurDate)
                    MW=self.GetMeanWeight(s,CurDate)
                    self.PredVal+=[Lambda('PredVal'+subscript,\
                        lambda  useArea=useArea,\
                                lnGrandMean=self.lnGrandMean,\
                                RBM=RBM,\
                                SE=SE,\
                                YE=YE,\
                                TE=TE,\
                                MW=MW:\
                                useArea*RBM/MW*exp(lnGrandMean+SE+YE+TE),trace=False )]
                    self.ObsVal+=[Poisson('ObsVal'+subscript,self.PredVal[-1],observed=True,value=CurNumAnimals)] 
            

if __name__ == "__main__":
  import pickle
  pfile="../MCMC/Jervis/Jervis.pickle"
  SiteNumber,Sites,VirginSites,CoastLength=pickle.load(open(pfile,"rb"))  
  
  ConstantValues.MinYear=1999  #Earliest year to be considered in calculations
  
  test=EFAProdModel(SiteNumber,Sites,VirginSites,CoastLength,oldB=1)
  #print(test.GetMeanWeight(16,[2005,2,14]))  

