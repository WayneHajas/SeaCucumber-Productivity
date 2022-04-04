'''
Class to represent harvest that has taken place at a particular site.
2018-02-14
  Harvest date had been forced to the first day of the month.  Got rid of that.

'''
from numpy import ndarray,ceil
import sys
sys.path.append('D:/Analyses/CukeEFAProd/pyfunctions')
from ADO import adoBaseClass as daoBaseClass
from ToDecimalYear import ToDecimalYear

class Harvest():
  def __init__(self,mdbfile,Survey,Site,MinYear=1995,MaxYear=sys.maxsize,ConvFactor=0.45359237):
    '''
    Harvest(mdbfile,Survey,Site,MinYear=1995,MaxYear=sys.maxsize,ConvFactor=0.45359237)
    Represents harvest that occurs from year MinYear to year MaxYear.
    
    The amount harvested is multiplied by ConvFactor.  So if the database gives 
    harvest in pounds, ConvFactor=0.45359237 will convert harvest to units of kilograms   
    
    '''
    self.Survey=Survey
    self.Site=Site
    self.MinYear=MinYear
    self.MaxYear=MaxYear
    self.ConvFactor=ConvFactor
    
    #Create query for retrieving the data
    query= 'select HarvestHeaders.Year,HarvestHeaders.FMonth ,HarvestHeaders.FDay ,  '
    query+='Sum(HarvestValidation.NetWeight)*  '
    query+=str(self.ConvFactor)
    query+=' AS HarvKg '
    query+='FROM (HarvestValidation INNER JOIN HarvestHeaders ON  '
    query+='HarvestValidation.HHFKey = HarvestHeaders.HHKey) INNER JOIN  '
    query+='HarvestProject ON HarvestHeaders.HPFKey = HarvestProject.HPkey '
    query+='GROUP BY HarvestProject.Project, HarvestHeaders.Year, HarvestHeaders.FMonth,HarvestHeaders.FDay,  '
    query+='HarvestHeaders.Site, HarvestHeaders.Year, HarvestHeaders.FMonth '
    query+='HAVING (((HarvestProject.Project)="'
    query+=self.Survey
    query+='") AND ((HarvestHeaders.Site)= '
    query+=str(self.Site)
    query+=')'
    query+=' and (HarvestHeaders.Year >= '+str(MinYear)+") "
    if MaxYear !=sys.maxsize:
        query+=' and (HarvestHeaders.Year <= '+str(MaxYear)+") "
    query+=') '
    query+='ORDER BY HarvestHeaders.Year,  '
    query+='HarvestHeaders.FMonth, '
    query+='HarvestHeaders.FDay;'

    dataSource=daoBaseClass(mdbfile,query)
    
    #Harvest dates are converted to decimal-years.
    self.hdata=[ [ToDecimalYear([t[0],t[1],t[2]]),t[3]]   for t in dataSource.GetALL()]
    del dataSource
    
    #The harvest dates as a separate list
    self.HarvestDates=[t[0] for t in self.hdata]
    
  def GetHarv(self,testdate):
      
    #Force date to be a decimal-year  
    DecDate=ToDecimalYear(testdate)
    
    #Deal with list of dates
    if isinstance(DecDate,(list,ndarray)):
        result=[self.HarvFromDecimalDate(t) for t in DecDate ]
        return(result)
    result=self.HarvFromDecimalDate(DecDate)
    return(result)
  
  def HarvFromDecimalDate(self,DecDate,fudge=.5/365.25):
    '''
    Harvest.HarvFromDecimalDate(DecDate,fudge=.5/365.25)
    Get harvest that occured at DecDate; +/-fudge
    DecDate is the date as a decimal year
    fudge is just the tolerance for specifying the date.
    '''
    if (DecDate<self.MinYear):return(0.)
    if (DecDate>ceil(self.MaxYear)):return(0.)
    
    for h in self.hdata:
      #Check if any of the harvest dates correspond to the testdate
      if (h[0]>=(DecDate-fudge)) and (h[0]<=(DecDate+fudge)):
          return(h[1])
          
    #No matches
    return(0.)

  def GetHarvestDate(self,MinYr=None,MaxYr=None):
      result=self.HarvestDates
      if MinYr:
          result=[t for t in result if t>=MinYr]
      if MaxYr:
          result=[t for t in result if t<=MaxYr]
      result.sort()
      return(result)
   
  		
    
if __name__ == "__main__":
  mdbfile='t:\SeaCuke_Bio.mdb'
  Survey='Laredo Inlet'
  Site=16
  Transect=135
  MinYear=1995
  ConvFactor=0.45359237
  from ToDecimalYear import FromDecimalYear
  
  h=Harvest(mdbfile,Survey,Site,MinYear=MinYear,ConvFactor=ConvFactor)
  print()
  for t in h.hdata:
    print (FromDecimalYear(t[0]),t[1])
  print()
 
  from ToDecimalYear import FromDecimalYear
  from datetime import datetime
  for t in h.hdata:
      print(FromDecimalYear(t[0]),t)
      
  print(datetime(2002,10,1), h.GetHarv(datetime(2002,10,1)))
  print(datetime(2002,11,1), h.GetHarv(datetime(2002,11,1)))

  