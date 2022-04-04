'''2019-03-08
        Virgin-transect chosen according to year.
        Assume there is only one survey during that year.
    2017-04-11
        Fix an error in DayTran.DepthProfInRange
        Get the last(shallowest) quadrat
    2017-08-04
        rewrote DayTran.DepthProfInRange to get the correct number of quadrats after trimming
    2017-08-23
        Create Transect.HaveEnoughCuke(LowBnd=5)
        Checks to see if LowBnd or more cukes were observed during any of the surveys
	2017-09-14
	
		modify DayTran.DepthProfInRange
		Resulting profile is from the first quadrat in the depth-range to the last quadrat depth-range.  
		Quadrats in the middle might not be in range.
		
'''

from operator import  attrgetter
from numpy import ndarray
import sys
from sys import maxsize as maxint
from ADO import adoBaseClass as daoBaseClass
from ToDecimalYear import ToDecimalYear,FromDecimalYear
from constants import ConstantValues

class DayTran():
  def __init__(self,mdbfile,Survey,Site,Transect,YearMonthDay,DateByMonth=False):
      
    #Force the date to a format
    dtDate= FromDecimalYear(YearMonthDay) 
    self.year,self.month,self.day=dtDate.year,dtDate.month,dtDate.day
    self.TranNum=Transect
    self.Site=Site
    query= 'SELECT cint((Densities.Distance)/5-1) AS Qnum, Densities.ChartDepth, '
    query+='Densities.CountLeft+Densities.CountRight AS NumCukes '
    query+='FROM Headers INNER JOIN Densities ON Headers.Key = Densities.HKey '
    query+='Where (((Headers.Project)="'
    query+=Survey
    query+='") AND ((Headers.EFA) Is Not Null) AND ((Densities.ChartDepth) Is Not Null) AND (Not((Densities.Distance) =0)) AND ((Headers.Site)= '
    query+=str(Site)
    query+=') AND ((Headers.Year)= '
    query+=str(self.year)
    query+=') AND  '
    query+='((Headers.Month)= '
    query+=str(self.month)
    query+=')AND '

    #ignore zero day-value
    query+='(Headers.Day= 0 or Headers.Day='
    query+=str(self.day)
    query+=')  AND '

    query+='((Headers.Transect)= '
    query+=str(Transect)
    query+=')) '
    query+='ORDER BY Densities.Distance;'

    dataSource=daoBaseClass(mdbfile,query)
    self.ydata=dataSource.GetALL()
    self.ydata.sort(key=lambda t:t[0]) #quadrats sorted according to distance-value in densities-table
    del dataSource
    
    nrec=len(self.ydata)
    for i in range(1,nrec):
      deltaD=self.ydata[i][0]-self.ydata[i-1][0]
      if deltaD>1:
        for j in range(self.ydata[i][0]+1,self.ydata[i-1][0]):
          frac=(j-self.ydata[i][0])/deltaD
          newQuad=[ j,\
                    self.ydata[i-1][1] +frac*(self.ydata[i][1] - self.ydata[i-1][1]),\
                    self.ydata[i-1][2] +frac*(self.ydata[i][2] - self.ydata[i-1][2])  ]
          self.ydata+=[newQuad]
    self.ydata.sort(key=lambda t:t[0]) #quadrats sorted according to distance quad number            
        
     
    self.GetDepthRange()
    
    if DateByMonth:
        if self.day>15:
            self.month+=1
            self.day=1
            if self.month>12:
                self.month=1
                self.year+=1
        else:
            self.day=1
    
   
  def GetDepthRange(self): 
  	self.DepthRange=[-maxint,maxint]
  	if len(self.ydata)>0:
  	  d=[t[1] for t in self.ydata]
  	  try:
  	    self.DepthRange=[min(d),max(d)]
  	  except:
  	    print('\ntransect 64')
  	    print(d)
  	    self.DepthRange=[min(d),max(d)]
  	
  def QuadInDepthRange(self,quad,DepthRange):
     try:
       if isinstance(quad[0],(list,ndarray)):
          ByQuad=[self.QuadInDepthRange(q,DepthRange) for q in quad]
          return(all(ByQuad))
     except:
       if isinstance(quad[0],(list,ndarray)):
          ByQuad=[self.QuadInDepthRange(q,DepthRange) for q in quad]
          return(all(ByQuad))
     result=(quad[1]>=DepthRange[0]) & (quad[1]<=DepthRange[1])
     return(result)
  		
      
  def DepthProfInRange(self,DepthRange):
    '''Generate the depth-profile from the first quadrat in the depth range to the last quadrat in the depth range.'''
    
    #Return empty list if there are no quadrats
    if not(self.ydata):
        return([])
    
    #Initial depth profile    
    depth=[t[1]  for t in self.ydata]
    
    #Take off initial quadrats that are out of depth-range
    while (depth[0]<DepthRange[0]) or (depth[0]>DepthRange[1]):
      depth.pop(0)
      
      #Expand depthrange to try to get at least one quadrat
      if not(depth):
          DepthRange2=[DepthRange[0]-0.25, DepthRange[1]+0.25   ]
          result=self.DepthProfInRange(DepthRange2)
          return(result)

    #Take off last quadrats that are out of depth-range
    while (depth[-1]<DepthRange[0]) or (depth[-1]>DepthRange[1]):
        depth.pop(-1)
      
        #Expand depthrange to try to get at least one quadrat
        if not(depth):
          DepthRange2=[DepthRange[0]-0.25, DepthRange[1]+0.25   ]
          result=self.DepthProfInRange(DepthRange2)
          return(result)
   #Quadrats in the middle that are out-of-range are left in the profile 
    
  
    return (depth)
    
  def MeasFit(self,Profile,index=None):
     '''Provide a measurement of how well the transect fits a profile.  A small value is better.
         index indicates the starting quadrat
         Default (index=None) is to generate a value for every quadrat in the transect'''
         
     if index==None:
        i=range(len(self.ydata)-len(Profile) +1  )
        result=[self.MeasFit(Profile,index=j) for j in i]
        return(result)
  	
     #invalid value for index    
     if (index>((len(self.ydata)-len(Profile)))):
         return(maxint)
        
     sq=[ (Profile[j]-self.ydata[index+j][1])**2      for j in range(len(Profile))]
     result=sum(sq)
     return(result)
  
  def GetMinSSQ(self,Profile):
      FitSSQ=self.MeasFit(Profile)
      result=min(FitSSQ)
      return(result)

  def GetAllProfileForLength(self,Length):
      nquad=len(self.ydata)
      result=[  [t[1] for t in self.ydata[i:i+Length]]   for i in range(nquad-Length+1)]
      return(result)

  
  def ReduceToFit(self,Profile):
    FitSSQ=self.MeasFit(Profile)
    try:
      minSSQ=min(FitSSQ)
      minIndex=FitSSQ.index(minSSQ)
    except:
      print('transect 116',Profile)
      minSSQ=min(FitSSQ)
      minIndex=FitSSQ.index(minSSQ)
    try:
      self.ydata=self.ydata[minIndex:minIndex+len(Profile)]
    except:
      print('transect 129',minIndex)
      self.ydata=self.ydata[minIndex:minIndex+len(Profile)]
    self.GetDepthRange()
    
  def Summary(self):
    nquad=len(self.ydata)
    ncuke=sum(list(map(lambda x:x[-1]  ,self.ydata)))
    result={'year':self.year,'month':self.month,'day':self.day,'DepthRange':self.DepthRange,'nquad':nquad,'ncuke':ncuke}
    return(result)
    
  def GetNumCuke(self):
  	s=self.Summary()['ncuke']
  	return(s)
   
  def GetArea(self):
      nquad=self.Summary()['nquad']
      global ConstantValues
      result=nquad*ConstantValues.QuadArea
      return(result)
  def GetTranLen(self):
      nquad=self.Summary()['nquad']
      global ConstantValues
      result=nquad*ConstantValues.QuadLen
      return(result)

  def TrimToDepthRange(self,DepthRange):
      
    n=len(self.ydata)
    TempCopy=[t +[self.QuadInDepthRange(t,DepthRange)] for t in self.ydata  ]
    partList=[]
    if TempCopy[0][-1]:
        partList=[[TempCopy[0]]]
    for i in range(1,n):
        if TempCopy[i-1][-1]  and  TempCopy[i][-1]:
            partList[-1]+=[TempCopy[i]]
        if not(TempCopy[i-1][-1])  and  TempCopy[i][-1]:
            partList+=[[TempCopy[i]]]
            
    partList.sort(key=lambda t:-len(t))
    self.ydata=[t[:-1]  for t in partList[0]] 
    self.GetDepthRange() 
    
  
    		
class Transect():    		
  def __init__(self,mdbfile,Survey,Site,TransectNumber,WithTrim=False,track=True,DateByMonth=False):
    if track:print(' TransectNumber ',Survey,Site,TransectNumber)
    self.mdbfile=mdbfile
    self.Survey=Survey
    self.Site=Site
    self.TransectNumber=TransectNumber
    self.GetDayOfSurvey()
    try:
        self.AllQuad=[DayTran(self.mdbfile,self.Survey,self.Site,self.TransectNumber,DOS,DateByMonth=DateByMonth) for DOS in self.DayOfSurvey]
   
    except:
        self.AllQuad=[DayTran(self.mdbfile,self.Survey,self.Site,self.TransectNumber,DOS,DateByMonth=DateByMonth) for DOS in self.DayOfSurvey]
    try:
      self.AllQuad=[aq for aq in self.AllQuad   if len(aq.ydata)>0]
    except:
      harvestObject.set_trace()
      self.AllQuad=[aq for aq in self.AllQuad   if len(aq.ydata)>0]
    
    #Reduce Day-of-Survey to remaining surveys    
    self.DayOfSurvey=[{'year':aq.year,'month':aq.month,'day':aq.day}  for aq in self.AllQuad]
    self.DayOfSurvey=sorted( self.DayOfSurvey,key=lambda k:(k['year'],k['month'],k['day']))
    self.DecimalYear=ToDecimalYear(self.DayOfSurvey)

    self.GetCommonDepthRange()
    try:
      if WithTrim:self.Trim()
    except:
      print ('transect 169')
      print ('self.TransectNumber',self.TransectNumber)
      if WithTrim:self.Trim()
      

  def GetDayOfSurvey(self):
    global ConstantValues
    query= 'SELECT Headers.Year, Headers.Month, Headers.Day '
    query+='FROM Headers INNER JOIN Densities ON Headers.Key = Densities.HKey '
    query+='GROUP BY Headers.Project, Headers.Site, Headers.Transect, Headers.Year, Headers.Month, Headers.Day '
    query+='HAVING ( '
    query+=     '(Headers.Project= "'
    query+=         self.Survey
    query+=         '") AND '
    query+=     '(Headers.Site= '
    query+=         str(self.Site)+' ) '
    query+=     ' AND (Headers.Transect= '
    query+=         str(self.TransectNumber)+')'
    query+=     ' AND (Headers.Year>= '
    query+=         str(ConstantValues.MinYear)+')'
    query+=     ' AND (Headers.Year<= '
    query+=         str(ConstantValues.MaxYear)+')'
    query+= ') '
    query+=' order by  Headers.Year,Headers.Month,Headers.Day   '    
    query+=';'
    
    try:
        dataSource=daoBaseClass(self.mdbfile,query)
    except:
        print ('149 ',query)
    self.DayOfSurvey=[ {'year':t[0],'month':t[1],'day':t[2]}  for t in dataSource.GetALL()]
    del dataSource

  def HasChartDepth(self,SurveyDate):
    #force date into datetime format.
    dtDate=FromDecimalYear(SurveyDate)
    year=dtDate.year
    month=dtDate.month
    day=dtDate.day
    query='SELECT Count(Densities.ChartDepth) AS NumChartDepth '
    query+='FROM Densities INNER JOIN Headers ON Densities.HKey = Headers.Key '
    query+='GROUP BY Headers.Project, Headers.Transect, Headers.Site, Headers.Year, Headers.Month '
    query+='Having (((Headers.Project)= "'
    query+=self.Survey
    query+='") AND ((Headers.Transect)= '
    query+=str(self.TransectNumber)
    query+=') AND ((Headers.Site)= '
    query+=str(self.Site)
    query+=') AND ((Headers.Year)= '
    query+=str(year)
    query+=') AND ((Headers.Month)= '
    query+=str(month)
    
    #ignore zero day-value
    query+=')  AND (Headers.Day=0 or Headers.Day= '
    query+=str(day)
    query+=') AND  '
    query+='(Count(Densities.ChartDepth) Is Not Null));'
    
    dataSource=daoBaseClass(self.mdbfile,query)
    NumQuad=dataSource.GetALL()
    del dataSource
    result=NumQuad[0][0]>0
    return(result)
  
 
  def GetCommonDepthRange(self):
  	MinDepth=[]
  	MaxDepth=[]
  	for y in self.AllQuad:
  	  d=y.DepthRange
  	  MinDepth+=[d[0]]
  	  MaxDepth+=[d[1]]
  	try:
  		self.CommonDepthRange=[max(MinDepth),min(MaxDepth)]
  	except:
  		self.CommonDepthRange=[-maxint,maxint]
  
  def DepthProfInRange(self):
      'Generate the depth-profile for each year and select the shortest one.'
      
      #Shortest transect-length associated with depth-range      
      dp=[ aq.DepthProfInRange(self.CommonDepthRange)   for aq in self.AllQuad]      
      dp=sorted(dp,key=lambda x:len(x))
      if dp==[]:
          return(None)
      targLen=len(dp[0])
      if targLen<=0:
          return([])
          
      #Check all the transect-segments of the target-length.  Get the best one
      result=self.GetBestProfile(targLen)
      return(result)

   
  def Trim(self):
  	try:
  		Profile=self.DepthProfInRange()
  	except:
  		print ('transect 292')
  		print ('self.TransectNumber',self.TransectNumber)
  		print ('self.Site',self.Site)
  		print ('self.CommonDepthRange',self.CommonDepthRange)
  		Profile=self.DepthProfInRange()
  	for y in self.AllQuad:
  		try:
  			y.ReduceToFit(Profile)
  		except:
  			print ('transect 282 ',self.TransectNumber)
  			y.ReduceToFit(Profile)
  			
  	if Profile !=None:
  	  self.nTrimQuad=len(Profile)
  	  self.GetCommonDepthRange()
  	else:
  	  self.nTrimQuad=0
  	  self.GetCommonDepthRange()
  	
  def Summary(self):
      result=[aq.Summary()  for aq in self.AllQuad]      
      result=sorted(result,key=lambda x:x['year'])
      return(result)
  
  def GetTranAtYear(self,fltD):
    '''Transect.GetTranAtYear(fltD)
       Give the transect at a specfic time-value.
       fltD is the survey-date; as a datetime,decimalyear or [y,m,d]'''
    dtDate=FromDecimalYear(fltD)
    
    #multiple dates    
    if isinstance(dtDate,list):
        result=[self.GetTranAtYear(dtd)   for dtd in dtDate]
        return(result)
    
    #Single date    
    y,m,d=dtDate.year,dtDate.month,dtDate.day
    result=[aq     for aq in self.AllQuad    if (aq.year==y) and (aq.month==m)and (aq.day==d) ]
    if len(result)!=1:
      print (' transect 300 ',fltD, self.TransectNumber)
      return(None)
    result=result[0]
    return(result)
    
  def GetNumCuke(self,fltD=None):
      TranAtYear=GetTranAtYear(fltD)
      
      #multiple surveys
      if isinstance(TranAtYear,list):
          result=[t.GetNumCuke()  for t in TranAtYear]
          return(result)
          
      #single survey
      result=TranAtYear.GetNumCuke()
      return(result)

  def HaveEnoughCuke(self,LowBnd=5):
      '''Determine if transect has at least LowBnd cukes in any of the surveys'''
      NumCuke=[t.GetNumCuke()  for t in self.AllQuad]
      NumCuke=[cuke for cuke in NumCuke if cuke!=None]
      maxCuke=max(NumCuke)
      result=(maxCuke>=LowBnd)
      return(result)
      

    

  def YearOfSurvey(self):
    result=list(map(lambda s:s['year'],self.DayOfSurvey))
    return(result)
    

  def GetMinSSQ(self,Profile):
      result=sum([t.GetMinSSQ(Profile)   for t in self.AllQuad])
      return(result)

  def GetAllProfileForLength(self,Length):
      result=[]
      for t in self.AllQuad:
          result+=t.GetAllProfileForLength(Length)
      return(result)
      
  def GetBestProfile(self,Length):
      AllProfile=self.GetAllProfileForLength(Length)
      SSQ=[ self.GetMinSSQ(t)   for t in AllProfile]
      index=SSQ.index(min(SSQ))
      result=AllProfile[index]
      return(result)
  
  def GetNumCuke(self,SurveyDate=None):
      '''transect.GetNumCuke(SurveyDate=None)
      if SurveyDate is un-specified, results for the first(virgin) survey are given'''
      
      if SurveyDate==None:
          useDate=FromDecimalYear(self.DayOfSurvey[0])
          result=self.GetNumCuke(useDate)
          return(result)
     
      useDate =FromDecimalYear(SurveyDate)
      if isinstance(useDate,(list,ndarray)):
          result=[self.GetNumCuke(ud)  for ud in useDate    ]
          return(result)
      dataDate=FromDecimalYear(self.DayOfSurvey)
      try:
         index=dataDate.index(useDate)
      except:
          #There was no survey on the date
          return(None)
      CurSurvey=self.AllQuad[index]
      result=CurSurvey.GetNumCuke()
      return(result)
   
  def GetArea(self,SurveyDate=None):
      '''transect.GetArea(SurveyDate=None)
      if SurveyDate is un-specified, results for the first(virgin) survey are given'''
      
      if SurveyDate==None:
          useDate=FromDecimalYear(self.DayOfSurvey[0])
          result=self.GetArea(SurveyDate=useDate)
          return(result)
     
      useDate =FromDecimalYear(SurveyDate)
      if isinstance(useDate,(list,ndarray)):
          result=[self.GetArea(ud)  for ud in useDate    ]
          return(result)
      dataDate=FromDecimalYear(self.DayOfSurvey)
      try:
         index=dataDate.index(useDate)
      except:
          #There was no survey on the date
          return(None)
      CurSurvey=self.AllQuad[index]
      result=CurSurvey.GetArea()
      return(result)
   
  def GetTranLen(self,SurveyDate=None):
      '''transect.GetTranLen(SurveyDate=None)
      if SurveyDate is un-specified, results for the first(virgin) survey are given'''
      
      if SurveyDate==None:
          useDate=FromDecimalYear(self.DayOfSurvey[0])
          result=self.GetTranLen(useDate)
          return(result)
     
      useDate =FromDecimalYear(SurveyDate)
      if isinstance(useDate,(list,ndarray)):
          result=[self.GetTranLen(ud)  for ud in useDate    ]
          return(result)
      dataDate=FromDecimalYear(self.DayOfSurvey)
      try:
         index=dataDate.index(useDate)
      except:
          #There was no survey on the date
          return(None)
      CurSurvey=self.AllQuad[index]
      result=CurSurvey.GetTranLen()
      return(result)

class VirginTransect(DayTran):
  def __init__(self,CurTran,vyear=None,DateByMonth=False):
      
    #Force the date to a format
    if vyear==None:
        dtDate= FromDecimalYear(CurTran.DayOfSurvey[0]).year
    else:
        dtDate= FromDecimalYear(vyear).year
    self.year   = dtDate
    self.month  = FromDecimalYear(CurTran.DayOfSurvey[0]).month
    self.day    = FromDecimalYear(CurTran.DayOfSurvey[0]).day
    
    if DateByMonth:
        if self.day>15:
            self.month+=1
            self.day=1
            if self.month>12:
                self.month=1
                self.year+=1
        else:
            self.day=1

    
    self.TranNum=CurTran.TransectNumber
    self.Site=CurTran.Site
    query= 'SELECT cint((Densities.Distance)/5-1) AS Qnum, Densities.ChartDepth, '
    query+='Densities.CountLeft+Densities.CountRight AS NumCukes '
    query+='FROM Headers INNER JOIN Densities ON Headers.Key = Densities.HKey '
    query+='Where (((Headers.Project)="'
    query+=CurTran.Survey
    query+='") AND ((Headers.EFA) Is Not Null) AND ((Densities.ChartDepth) Is Not Null) AND (Not((Densities.Distance) =0)) AND ((Headers.Site)= '
    query+=str(CurTran.Site)
    query+=') AND ((Headers.Year)= '
    query+=str(self.year)
    query+=') AND  '

    query+='((Headers.Transect)= '
    query+=str(self.TranNum)
    query+=')) '
    query+='ORDER BY Densities.Distance;'

    dataSource=daoBaseClass(CurTran.mdbfile,query)
    self.ydata=dataSource.GetALL()
    self.ydata.sort(key=lambda t:t[0]) #quadrats sorted according to distance-value in densities-table
    del dataSource
    
    nrec=len(self.ydata)
    for i in range(1,nrec):
      deltaD=self.ydata[i][0]-self.ydata[i-1][0]
      if deltaD>1:
        for j in range(self.ydata[i][0]+1,self.ydata[i-1][0]):
          frac=(j-self.ydata[i][0])/deltaD
          newQuad=[ j,\
                    self.ydata[i-1][1] +frac*(self.ydata[i][1] - self.ydata[i-1][1]),\
                    self.ydata[i-1][2] +frac*(self.ydata[i][2] - self.ydata[i-1][2])  ]
          self.ydata+=[newQuad]
    self.ydata.sort(key=lambda t:t[0]) #quadrats sorted according to distance quad number            
        
     
    self.GetDepthRange()  

if __name__ == "__main__":
  mdbfile='t:\SeaCuke_Bio.mdb'
  Survey='Laredo Inlet'
  Site=8
  TransectNumber=5
  YearMonth=[1999,2,18]
  WithTrim=True
  
  global ConstantValues
  ConstantValues.MinYear=1998

  
  
  t=Transect(mdbfile,Survey,Site,TransectNumber,WithTrim=WithTrim)
  print('transect 284')
  Profile=t.DepthProfInRange()
  print('transect 286')
  v=VirginTransect(t)
  print(v.Summary())

  print()
  print(t.GetBestProfile(3))
  print()
  for s in t.AllQuad:
      d=[r[1] for r in s.ydata]
      print(d)
  print()  
  for s in t.AllQuad:
      print(s.Summary()['ncuke'])
  
  
  print ('done')
  
  

