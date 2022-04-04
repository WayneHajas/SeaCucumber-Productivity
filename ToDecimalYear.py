from datetime import datetime as dt
from numpy import ndarray
from calendar import monthrange

#An acknowledged uncertainty when dates are converted to decimals and back
epsday=0.5/365.25

def ToDecimalYear(dtdate):
    if isinstance(dtdate,(list,ndarray)):
        if (len(dtdate)==3) and all ([isinstance(t,int) for t in dtdate ]):
            #Force day-of-month to be at least one
            fixday=max([1,dtdate[2]])
            result=ToDecimalYear(dt(dtdate[0],dtdate[1],fixday))
            return(result)
        result=[ToDecimalYear(t) for t in dtdate]
        return(result)
    
    if isinstance(dtdate,dict):
        result=ToDecimalYear([ dtdate['year'],dtdate['month'],dtdate['day']  ])
        return(result)
    
    #A single value
    if isinstance(dtdate,(float,int)):
        return(dtdate)
    
    try:
        y=dtdate.year
    except:
        print('ToDecimalYear 29',dtdate)
        y=dtdate.year
    d    =dtdate.toordinal()
    dprev=dt(y,1,1).toordinal()
    dafte=dt(y+1,1,1).toordinal()
    result=y+(d-dprev)/(dafte-dprev)
    return(result)

def FromDecimalYear(fltD):
    if isinstance(fltD,dt):
        return(fltD)
    if isinstance(fltD,(list,ndarray)):
        if len(fltD)==3:
            
            #Date is given as three integers
            if isinstance(fltD[0],int) and isinstance(fltD[1],int) and isinstance(fltD[2],int): 
                try:
                    result=dt( fltD[0],max(fltD[1],1),max([fltD[2],1])  )
                except:
                    print('\n ToDecimalYear 39')
                    print(fltD)
                    result=dt( fltD[0],fltD[1],fltD[2]  )
                return(result)
        result=[FromDecimalYear(t) for t in fltD]
        return(result)
    
    #Date is given as an integer-year    
    if(isinstance(fltD,int)):
        try:
            return(dt(fltD,1,1))
        except:
            print('\nToDecimalYear 49')
            print(fltD)
            return(dt(fltD,1,1))
        
    #Date occurs as a year-month-day dictionary
    if isinstance(fltD,dict):
        #force day of month to be at least one
        fixday=max([1,fltD['day']])
        fixmonth=max([1,fltD['month']])
        try:
          return(dt(fltD['year'], fixmonth, fixday ))
        except:
          return(dt(fltD['year'], fltD['month'], fixday ))
    
    #Date is a whole-year    
    if fltD==int(fltD):
        return(dt(int(fltD),1,1))
    
    y=int(fltD)
    ByMonth=[ToDecimalYear(dt(y,i,1)) for i in range(1,13)]
    m=sum(t <=fltD for t in ByMonth)
    #ByDay=[epsday+ToDecimalYear(dt(y,m,i)) for i in range(1,1+monthrange(y,m)[-1])   ]
    #day=sum(t <=fltD for t in ByDay)
    ssq=[ [i,(fltD-ToDecimalYear(dt(y,m,i)))**2] for i in range(1,1+monthrange(y,m)[-1])   ]
    minssq=min(ssq,key=lambda t:t[1])
    day=minssq[0]
    
    
    
    try:
        return(dt(y,m,day))
    except:
        print('ToDecimalYear 42')
        print(y,m,day,fltD)
        return(dt(y,m,day))

def CombineDateLists(dL1, dL2):
    '''CombineDateLists(dL1, dL2)
    Assume dL1 and dL2 are both lists'''
    catlist=ToDecimalYear(dL1)+ToDecimalYear(dL2)
    
    #get rid of duplicates
    trunclist=list(set(catlist))
    trunclist.sort()
    return(trunclist)
    
def GetCalendarYear(DateValue):
    if isinstance(DateValue,(list,ndarray)):
        if len(DateValue)==3:
            if isinstance(DateValue[0],int) and isinstance(DateValue[2],int) and isinstance(DateValue[2],int): 
                return(DateValue[0])
        result=[GetCalendarYear(t) for t in DateValue]
        return(result)
        
    if isinstance(DateValue,float):
        return(int(DateValue))
    if isinstance(DateValue,int):
        return(DateValue)
    if isinstance(DateValue,dict):
        return(DateValue['year'])
    if isinstance(DateValue,dt):
        return(DateValue.year)

if __name__ == "__main__":
    dtdate=[dt(2005,12,31),dt(2005,12,1),dt(2005,1,31),dt(2005,1,1),dt(2005,9,17)]
    TestResult=ToDecimalYear(dtdate)
    TestResult2=FromDecimalYear(TestResult)
    
    for i in range(len(dtdate)):
        print(dtdate[i],TestResult[i],TestResult2[i])
