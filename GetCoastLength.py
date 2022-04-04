from numpy import ndarray
from ADO import adoBaseClass as daoBaseClass

def GetCoastLength(Survey,Site=None,mdbfile='S:\\analyses\\EFA_Productivity_Model.20180929\\data\\CoastLength.accdb'):
    if Site==None:
        result=GetCoastLength(Survey,Site=[0,2,4,8,16],mdbfile=mdbfile)
        return(result)
    if isinstance(Site,(list,ndarray)):
        result={}
        for s in Site:
            result[s]=GetCoastLength(Survey,Site=s,mdbfile=mdbfile)
        return(result)

    #Get Coastlength for combination of Survey and site
    query= 'SELECT CoastLength.CoastLength_metres '
    query+='FROM CoastLength '
    query+='WHERE (((CoastLength.Location)="'
    query+=Survey
    query+='") AND ((CoastLength.Site)= '
    query+=str(Site)
    query+='));'
    dataSource=daoBaseClass(mdbfile,query)
    result=dataSource.GetVariable('CoastLength_metres')[0]
    return(result)
  
if __name__ == "__main__":
    t1=GetCoastLength('Jervis Inlet',Site=0)
    print(t1)
    t2=GetCoastLength('Jervis Inlet',Site=[0,2])
    print(t2)
    t3=GetCoastLength('Jervis Inlet')
    print(t3)