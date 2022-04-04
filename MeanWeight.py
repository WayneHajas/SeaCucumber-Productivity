

'''20160307
    Include 95%CB, +/-2StErr, as a statistic for estimated mean weight.#run D:/Analyses/CukeEFAProd/pyfunctions/MeanWeight
2016-04-13
    From EFA or non-EFA, group together all animals from the same calendar year. 
    
    Keep the lowest month and lowest day-values to generate a date.  For now, this will work because 
    data is always collected during the same month.  Could cause problems later.

20210125
    For some reason, ACCESS no longer likes the CDbl operator.  I am just taking it out.

'''

from numpy import ndarray,floor,sqrt,empty


from ADO import adoBaseClass as daoBaseClass
from ToDecimalYear import ToDecimalYear,FromDecimalYear


class MeanWeight():
  def __init__(self,mdbfile,Survey,Site=None):
     '''MeanWeight(mdbfile,Survey,Site)
     Source of mean-weight values for a site.
     Node_MeanWeight.GetMeanWeight(Date) will get closest(in time) mean weight'''  

     self.mdbfile=mdbfile
     self.Survey=Survey
     self.Site=Site
     
     #Average Weights from EFAbiosamples   
     self.ReadEFAAvg()
     
     #Non EFAsamples     
     self.ReadNonEFAAvg()

       
     
  def ReadEFAAvg(self):
      #Note that weight in grams is convert to kilograms          
      query= 'SELECT  HarvestHeaders.Year, min(HarvestHeaders.FMonth) as FMonth, min(HarvestHeaders.FDay) as FDay, '
      query+='Avg(EFABioSample.DrainedWgt)*0.001 AS AvgWgt,    '
      query+='StDev((EFABioSample.DrainedWgt))/Sqr(Count(EFABioSample.DrainedWgt))*0.001 AS StErrWgt '
      
      query+='FROM ((HarvestDives INNER JOIN EFABioSample ON HarvestDives.HDKey = EFABioSample.HDKey)  '
      query+='INNER JOIN HarvestHeaders ON HarvestDives.HHFKey = HarvestHeaders.HHKey)  '
      query+='INNER JOIN HarvestProject ON HarvestHeaders.HPFKey = HarvestProject.HPkey  '
      
      query+='GROUP BY HarvestHeaders.Year,  HarvestProject.Project, HarvestHeaders.Site '
      query+='HAVING ((HarvestProject.Project="'
      query+=self.Survey+'")'
      if self.Site!=None:
        query+='AND (HarvestHeaders.Site='
        query+=str(self.Site)+') '
      query+=' AND ((Avg(EFABioSample.DrainedWgt))<1000)  '
      query+=') '
      query+='order by HarvestHeaders.Year ;'
           
      dataSource=daoBaseClass(self.mdbfile,query)
      EFAavgWgt=dataSource.GetALL()
      self.EFAavgWgt=[{'year':t[0],'month':t[1],'day':t[2],'date':FromDecimalYear(t[:3]),\
                  'AvgWgt':t[3],'StErr':t[4],'tau':1/t[4]/t[4],'95CB':[t[3]-2*t[4],t[3]+2*t[4]]}   for t in EFAavgWgt]
      self.nEFA=len(self.EFAavgWgt)

  def ClosestEFA(self,year) :   
       decimalyear=ToDecimalYear(year)
                  
       SSQ=[ (decimalyear-ToDecimalYear(y['date']))**2  for y in self.EFAavgWgt]
       index=SSQ.index(min(SSQ))
       return(self.EFAavgWgt[index])
        
  def ReadNonEFAAvg(self):
      #Note that weight in grams is convert to kilograms
      query= 'SELECT Headers.Year, min(Headers.Month) as mmm,min(Headers.Day) as ddd, '
      query+='Avg(BioSample.DrainedWgt)*.001 AS AvgWgt,  '
      query+='StDev((BioSample.DrainedWgt))/Sqr(Count(BioSample.DrainedWgt))*.001 AS StErrWgt '
      query+='FROM BioSample INNER JOIN Headers ON BioSample.HKey = Headers.Key '
      query+='GROUP BY Headers.Project, Headers.Site, Headers.Year  '
      query+='HAVING ((Headers.Project="'
      query+=self.Survey+'") '
      if self.Site!=None:
        query+=" AND (Headers.Site= "
        query+=str(self.Site)+')'
      query+=') '
      query+='ORDER BY Headers.Year;'
    
      try:
          dataSource=daoBaseClass(self.mdbfile,query)
          NonEFAavgWgt=dataSource.GetALL()
          NonEFAavgWgt=[t for t in  NonEFAavgWgt if  not(any([s is None for s in t]))  ]
          self.NonEFAavgWgt=[{'year':t[0],'month':t[1],'day':t[2],'date':FromDecimalYear(t[:3]),\
                      'AvgWgt':t[3],'StErr':t[4],'tau':1/t[4]/t[4],'95CB':[t[3]-2*t[4],t[3]+2*t[4]]}   for t in NonEFAavgWgt]
          self.nEFA=len(self.NonEFAavgWgt)
      except:
          print('MeanWeight 93')
          print(query)
          dataSource=daoBaseClass(self.mdbfile,query)
          NonEFAavgWgt=dataSource.GetALL()
          NonEFAavgWgt=[t for t in  NonEFAavgWgt if  not(any([s is None for s in t]))  ]
          print(NonEFAavgWgt)
          self.NonEFAavgWgt=[{'year':t[0],'month':t[1],'day':t[2],'date':FromDecimalYear(t[:3]),\
                      'AvgWgt':t[3],'StErr':t[4],'tau':1/t[4]/t[4],'95CB':[t[3]-2*t[4],t[3]+2*t[4]]}   for t in NonEFAavgWgt]
          self.nEFA=len(self.NonEFAavgWgt)
          
      

  def ClosestNonEFA(self,year):
       decimalyear=ToDecimalYear(year)
                  
       SSQ=[ (decimalyear-ToDecimalYear(y['date']))**2  for y in self.NonEFAavgWgt]
       index=SSQ.index(min(SSQ))
       return(self.NonEFAavgWgt[index])
       
     

  def GetAvgWeight(self,year) :
       #Prefer to use EFA data
       decimalyear=ToDecimalYear(year)
       
       #Check for multiple years
       if isinstance(decimalyear,(list,ndarray)):
           result=[self.GetAvgWeight(t)  for t in decimalyear  ]
           return(result)
       
       dtEFA=self.ClosestEFA(year)
       dtEFA['SSQ']=(decimalyear-ToDecimalYear(dtEFA['date']))**2
       #An EFA-average is good enough if it occurs within 3 months.
       if dtEFA['SSQ']<(1./4./4.):
           return (dtEFA)
           
       #check the best non-EFA    
       dtnonEFA=self.ClosestNonEFA(year)
       dtnonEFA['SSQ']=(decimalyear-ToDecimalYear(dtnonEFA['date']))**2
       
       #Not a good average from either source
       if (dtnonEFA['SSQ']>1) and (dtEFA['SSQ']>1):
           print('\n MeanWeight 135')
           print('Did not get a good average weight for ',FromDecimalYear(year),' - site ',self.Site)
           if dtnonEFA['SSQ']<dtEFA['SSQ']:
               print('using ',dtnonEFA)
           else:
               print('using ',dtEFA)
               
       if dtnonEFA['SSQ']<dtEFA['SSQ']:
           return(dtnonEFA)
       
       return (dtEFA)   
       
  

if __name__ == "__main__":
  mdbfile='D:\Analyses\CukeNonParamProd\SeaCuke_Bio.mdb'
  Survey='Jervis Inlet'
  Site=0
  Year=1998
  Month=9
  day=15
  test=MeanWeight(mdbfile,Survey,Site)

  print('\n test.EFAavgWgt \n',test.EFAavgWgt)
  print('\n test.NonEFAavgWgt \n',test.NonEFAavgWgt)
  print('\n test.ClosestEFA(1999) \n',test.ClosestEFA(1999))
  print('\n test.ClosestNonEFA(1999) \n',test.ClosestNonEFA(1999))
  print('\n test.GetAvgWeight(1999) \n',test.GetAvgWeight(1999))
  print('\n test.GetAvgWeight(2000) \n',test.GetAvgWeight(2000))
  print('\n test.GetAvgWeight(2001) \n',test.GetAvgWeight(2001))
  print('\n test.GetAvgWeight([1962,8,5]) \n',test.GetAvgWeight([1962,8,5]))
  
