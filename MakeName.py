'''
MakeName.py
2017-01-24
'''
import datetime
from numpy import ndarray


def MakeName(prefix, lAttributes):
    '''
    MakeName(prefix, lAttributes)
    Make a suitable node-name.  Starts with prefix.  lAttributes-values are seperated but underscore'''
    
    result=prefix
    
    if not(isinstance(lAttributes,(list,ndarray))):
        result+='_'+sMakeName(lAttributes)
        return(result)
    
    for t in lAttributes:
        result+='_'+sMakeName(t)
    return(result)
 
def sMakeName( lAttributes):
   
    if isinstance(lAttributes,str):
        result=lAttributes
        return(result)
    
    if isinstance(lAttributes,int):
        result=str(lAttributes)
        return(result)
    
    if isinstance(lAttributes,datetime.datetime):
        result=str(lAttributes.year)+'_'+str(lAttributes.month)+'_'+str(lAttributes.day)
        return(result)
     
    #Just assume it can be converted to a suitable string
    return(str(lAttributes))
    
        
if __name__ == "__main__":
    prefix='xxx'
    
    x1=1
    x2='1'
    x3=datetime.datetime(2016,1,24)
    print(MakeName(prefix,x1))
    print(MakeName(prefix,x2))
    print(MakeName(prefix,x3))
    print()
    print(MakeName(prefix,[x1,x3]))
    print(MakeName(prefix,[x2,x3]))
    print(MakeName(prefix,[x1,x2,x3]))

    
    