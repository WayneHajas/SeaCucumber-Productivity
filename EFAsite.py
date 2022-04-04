'''2016-03-08
    Use year-of-first-survey for the transect in order to select virgin-transects for the site
    
2017-02-06
    In GetEventTimes, give the user the option of specifying a maximum  time-increment.
    For example:
        One event time is 1.7 and the next is 2.15 and the maximum time-increment is 0.1
            Would require 5 time-increments at 0.1
            Use 0.09 as time-increment.  Still require 5 increments.
2017-03-13
    Add a function, GetSurveyYear, which gives the calendar-years when surveys occur.
    
2007-04-04
    Modify GetEventTimes
    Only consider harvests if they occur after the last survey
2017-08-23
    Create Site.RemoveLowPopTransect(LowBnd=5)
    Remove all transects where the maximum number of cukes is less than LowBnd

2018-01-23
    add two more methods
    
    GetNumQuads gives the number of number of quadrats in the trimmed quadrats.  In a form compatible with the mimicWinBugs version of the model
    GetNumAnimals gives the number of animals in the trimmed quadrats.  In a form compatible with the mimicWinBugs version of the model
    '''

import sys
from numpy import ndarray
import numpy
from ADO import adoBaseClass as daoBaseClass


from transect import Transect, VirginTransect
from Harvest import Harvest
from constants import ConstantValues
from UniqueDates import UniqueDates
from ToDecimalYear import ToDecimalYear,FromDecimalYear,CombineDateLists
from MeanWeight import MeanWeight


class Site():
  def __init__(self,mdbfile,Survey,SiteNumber,CoastLength,DateByMonth=False):
    '''
    Site(mdbfile,Survey,SiteNumber,CoastLength,DateByMonth=False)
    
    Class to represent an EFA-site
    If DateByMonth, survey-dates will be shifted to the nearest first-day of the month.  That will reduce the number of Survey-dates.
    
    '''
    print('\n ')
    self.mdbfile=mdbfile
    self.Survey=Survey
    self.SiteNumber=SiteNumber
    self.CoastLength=CoastLength
    
    self.SiteHarv=Harvest(self.mdbfile,self.Survey,self.SiteNumber,MinYear=ConstantValues.MinYear,MaxYear=ConstantValues.MaxYear)
    self.ReadTransect(DateByMonth=DateByMonth)
    self.MeanWeight=MeanWeight(mdbfile,Survey,Site=self.SiteNumber)
    self.DateByMonth=DateByMonth
    
  def ReadTransect(self,track=True,DateByMonth=False):
    global ConstantValues
    query='SELECT DISTINCT Headers.Transect '
    query+='FROM Densities INNER JOIN Headers ON Densities.HKey = Headers.Key '
    query+='WHERE (((Headers.Project)="'
    query+=self.Survey
    query+='") AND ((Headers.EFA) Is Not Null) '
    query+=  ' AND (Headers.Year>= '
    query+=         str(ConstantValues.MinYear)+')'
    query+=  ' AND (Headers.Year<= '
    query+=         str(ConstantValues.MaxYear)+')'
    query+=') AND ((Headers.Site)= '
    query+=str(self.SiteNumber) 
    query+=') '
    query+='order by Headers.Transect ;'

    dataSource=daoBaseClass(self.mdbfile,query)
    
    #Transect numbers
    self.TranNum=dataSource.GetVariable('Transect')
    del dataSource
    
    #Create transects.  Each transect will be truncated to have the same lenght and similar depth-profile every time it is surveyed.    
    self.SurvTran=list(map(lambda t: Transect(self.mdbfile,self.Survey,self.SiteNumber,t,WithTrim=True,track=track,DateByMonth=DateByMonth),  \
                      self.TranNum))
    self.SurvTran=list(filter(lambda t:len(t.AllQuad)>0,self.SurvTran))
    
    #Dates when surveys occur.
    self.SetDayOfSurvey()
    
    #Create virgin-transect corresponding to the first Survey-date.  No truncation
    self.VirgTran=[ VirginTransect(t,vyear=self.VirginYear,DateByMonth=DateByMonth)  for t in self.SurvTran]    
    self.VirgTran=[t  for t in self.VirgTran if len(t.ydata)>0]
    
  
  def SetDayOfSurvey(self):
  	self.DayOfSurvey=[]
  	for t in self.SurvTran:
          self.DayOfSurvey+=t.DayOfSurvey
  	self.DayOfSurvey=UniqueDates(self.DayOfSurvey)
      
  	self.nSurvey=len(self.DayOfSurvey)
  	self.YearOfSurvey=list(map(lambda dos:  ToDecimalYear(dos)  ,self.DayOfSurvey))
  	self.VirginYear=int(self.YearOfSurvey[0])

     
  def GetTransect(self,tn,SurveyDate=None):
    '''
    Site.GetTransect(tn,SurveyDate=None)
    Retrieve a transect.
    tn is the transect-number
    If SurveyDate is defined, then result is an instance of DayTran
    If SurveyDate is undefined, result is an instance of Transect
    '''
    CurTran=[t for t in self.SurvTran if   t.TransectNumber==tn   ][0]
    if SurveyDate==None:
        return(CurTran)
    else:
        result=CurTran.GetTranAtYear(SurveyDate)
        return(result)

  def GetHarv(self,testDate):
      '''
      Site.GetHarv(testDate)
      Get harvest that occured on a specific date.
      '''
      result=self.SiteHarv.GetHarv(testDate)
      return(result)
        
  def GetEventTimes(self, MaxTimeInc=sys.maxsize):
      '''
      Site.GetEventTimes(MaxTimeInc=sys.maxsize)
      Get time values where survey and/or harvest occurs.

      If there is an increment greater than  MaxTimeInc between events, more time-values
      will be included so the maximum time-increment is never greater than MaxTimeInc    
      '''
      
      #Time Values that occur in data
      MaxTime=max(self.DayOfSurvey)#Date of last survey
      raw=[t for t in self.SiteHarv.HarvestDates if t<=MaxTime ]#Harvests before last survey
      raw+=self.DayOfSurvey#Combine harvest and survey dates.  Reduce to sorted unique values
      raw=list(set(raw))
      raw.sort()
      
      #Use these values if the maximum time-increment is set to its maximum
      if MaxTimeInc>=sys.maxsize:
          return(raw)
      if not(raw):
          return(raw)
      if len(raw) ==1:
        return(raw)
      
      #Build a list of event-times.  Where required, add more time-values to 
      #satisfy MaxTimeInc
      oldt=raw[0]
      result=[raw[0]]
      for newt in raw[1:]:
          #Do not need to break down the time increment
          if (newt-oldt)<=MaxTimeInc:
              result+=[newt]
          else:
              nincr=int(numpy.ceil((newt-oldt)/MaxTimeInc))
              delta=(newt-oldt)/nincr
              result+=[ oldt+delta*i for i in range(1, nincr+1)]
          oldt=newt
      return(result)
      
  def GetSurveyYears(self):
      '''Get calendar years, as integers where surveys took place'''
      #Force day-of-survey to be a datetime object
      x=FromDecimalYear(self.DayOfSurvey)
      
      y=[t.year for t in x]
      result=list(set(y))
      result.sort()
      return(result)
  
  def GetBestAvgWeight(self,year):
      '''result is a list of the form:=[{'year':'month':'day':,'date':FromDecimalYear(),\
                      'AvgWgt':,'StErr':,'tau':,'95CB':[]}
      '''
      result=self.MeanWeight.GetAvgWeight(year)
      return(result)
 
  def RemoveLowPopTransect(self,LowBnd=5):
        self.SurvTran=[t for t in self.SurvTran if t.HaveEnoughCuke(LowBnd=LowBnd)] 
        self.TranNum=[t.TransectNumber  for t in self.SurvTran]
    
        #Dates when surveys occur.
        self.SetDayOfSurvey()
    
        #Create virgin-transect corresponding to the first Survey-date.  No truncation
        self.VirgTran=[ VirginTransect(t,vyear=self.VirginYear,DateByMonth=self.DateByMonth)  for t in self.SurvTran]    
        self.VirgTran=[t  for t in self.VirgTran if len(t.ydata)>0]
      
  def GetNumQuads(self):     
      '''
      Site.GetNumQuads()
      Gives the number of quadrats in the trimmed version of each transect.  In a form that is compatible with EFAProdModel_MimicWinBugs
      '''
      result=[ t.nTrimQuad  for t in self.SurvTran]
      return(result)
  def GetNumAnimals(self,y=None):    
      SiteYears=y
      if not(SiteYears):
          SiteYears=self.GetSurveyYears()
      SurveyDates=self.DayOfSurvey
      
      result=[]
      for t in self.SurvTran:
          CurTran=[None for y in SiteYears]
          WithData=t.GetTranAtYear(SurveyDates)
          WithData=[s for s in WithData if s!=None]
          
          for w in WithData:
              y=w.year
              i=SiteYears.index(y)
              ncuke=w.GetNumCuke()
              CurTran[i]=ncuke
          result+=[CurTran]
      return(result)    
                  
      
        
class VirginSite(Site):
    def __init__(self,oriSite):
        '''VirginSite(oriSite)
        A classe to represent a site in its virgin-site.
        oriSite is an instantiation of the Site class.'''
        
        vDate=oriSite.DayOfSurvey[0]
        self.mdbfile=oriSite.mdbfile
        self.Survey=oriSite.Survey
        self.SiteNumber=oriSite.SiteNumber
        self.CoastLength=oriSite.CoastLength
        self.DayOfSurvey=[oriSite.DayOfSurvey[0]]
        
        self.SurvTran=oriSite.VirgTran
        self.TranLenStat()
        
    def TranLenStat(self):
        TL=[t.GetTranLen( ) for t in self.SurvTran]
        TL=[t for t in TL if t!=None]
        self.AvgTranLen=numpy.average(TL)
        stdTranLen=numpy.std(TL)
        n=len(TL)
        self.sterrTranLen=stdTranLen/numpy.sqrt(n)

if __name__ == "__main__":
  mdbfile='D:\Analyses\CukeNonParamProd\SeaCuke_Bio.mdb'
  Survey='Jervis Inlet'
  SiteNum=8

  global ConstantValues
  ConstantValues.MinYear=1995
  
  CoastLength=1000
  s=Site(mdbfile,Survey,SiteNum,CoastLength)
  v=VirginSite(s)

  y=s.GetSurveyYears()
 
  for t in s.SurvTran:
      year=[ty.Summary()['year']  for ty in t.AllQuad]
      ncuke=[ty.Summary()['ncuke']  for ty in t.AllQuad]
      nquad=[ty.Summary()['nquad']  for ty in t.AllQuad]
      print(t.TransectNumber,year,nquad[0],ncuke)
  print()
  s.RemoveLowPopTransect(LowBnd=5)
  for t in s.SurvTran:
      year=[ty.Summary()['year']  for ty in t.AllQuad]
      ncuke=[ty.Summary()['ncuke']  for ty in t.AllQuad]
      nquad=[ty.Summary()['nquad']  for ty in t.AllQuad]
      print(t.TransectNumber,year,nquad[0],ncuke)
  
  
  print ('done')

