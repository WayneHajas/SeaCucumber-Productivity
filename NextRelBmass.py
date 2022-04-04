'''
2017-01-23
NextRelBmass
Library to calculate next relative-biomass as per WCH productivity function
Harvest, if it occurs, is assumed to occur at the very beginning of the time-increment. Units are the same as for virgin biomass

2017-02-08
Add a wrapper so that calculations always start from a presumed virgin-state - prior to the first harvest.  
This is less computationally efficent than going from one reference time-point to the next, but it does allow for a simple structure when these calculations are incorporated into a Bayesian model.

2017-03-20
Corrected a serious bug in CurRelBmass
Previously, intitial year was based upon the first amount of harvest.  Now it is based upon the date of first harvest.
OldTime=useHarvestHist[0][0] is correct

2017-03-24
Corrected an error in ApplyBounds.
Now upper bound is applied.

2017-03-27
    Another limit on the maximum size of the time-increment.
    Corresponding biomass-increment must be less than 1/4 of B
    Corresponding biomass-increment must be less than 1/4 of 1-B
    
2017-07-26
    Restructure calculation of derivatives
2018-07-06
    Split the time-increment in two if a negative change is predicted for biomass    
    Split the time-increment in two if it is greater than one year
 2018-12-13
    Modify the calculation of relative biomass after harvest so that minimum possible vaue is  minRelBmass  
2019-04-04
    Modify the calculation of relative biomass after harvest 
        If post-harvest biomass is less than minRelBmass, don't calculate productivity
        If post-harvest biomass is less than 0.001, set value to 0.001
'''
from numpy import ndarray,abs
import warnings

from Harvest import Harvest
from ToDecimalYear import ToDecimalYear

#Target maximum allowable error for any time-increment when the WCH productivity model is applied.
#If the estimated error is larger than this value, the time-increment is split in half.
AllowableIntegrationError=.01

#If the time-increment is less than this value (in years), it is assumed there is no change in biomass.
MinTimeInc=.0001 #If time increment is smaller than this, produtivity is assume to be zero


def NextRelBmass(OldBiomass=None,OldTime=None,NewTime=None,a=None,b=None,xmax=None,fmax=None,VBiomass=None,CurHarvest=None,minRelBmass=1.e-3,maxRelBmass=1-1e-3):
    '''
    NextRelBmass(OldBiomass=None,OldTime=None,NewTime=None,a=None,b=None,xmax=None,fmax=None,VBiomass=None,CurHarvest=None,minRelBmass=1.e-3,maxRelBmass=1-1e-3)
    Calculate relative-biomass at time=NewTime as per WCH productivityh function
    Harvest, if it occurs, is assumed to occur at the very beginning of the time-increment
    
    * OldBiomass is the relative-biomass at the beginnin of the time-interval
    * OldTime is the decimal year at the beginning of the time-interval
    * NewTime is the decimal year at the end       of the time-interval
    * a,b,xmax,fmax are the parameters of the WCH productivity model 
    * VBiomass is the virgin biomass
    * CurHarvest is the amount of harvest  to occur at the beginning of the time-interval. Same units as VBiomass
    * minRelBmass and maxRelBmass are the minimum and maximum values of relative-biomass that will be considered.

    '''
    
    #There is an array of new-times.  Only allowe one harvest and it is at the very beginning of the first time-increment.
    if isinstance(NewTime,(list,ndarray)):
        result=[NextRelBmass(OldBiomass=OldBiomass,OldTime=OldTime,NewTime=NewTime[0],a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=CurHarvest,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)]
        ot=NewTime[0]
        for nt in NewTime[1:]:
            result+=[NextRelBmass(OldBiomass=result[-1],OldTime=ot,NewTime=nt,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=0,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)]
            ot=nt
        return(result)
        
    #A single new time
    
    #If OldBiomass is a node, use the value
    fOldBiomass=OldBiomass
    if not(isinstance(fOldBiomass,(float,int))):
        fOldBiomass=OldBiomass.value
    
    #There is a harvest
    if CurHarvest and (CurHarvest>0):
        ob=fOldBiomass-CurHarvest/VBiomass#Convert harvest to fraction of virgin-biomass and calculate relative-biomass after harvest
        if ob<=minRelBmass:
            result=max([ob,0.001])
            return(result)
        
        #Apply the time-increment starting from the post-harvest biomass        
        result=NextRelBmass(OldBiomass=ob,OldTime=OldTime,NewTime=NewTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=None,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
        return(result)
     
    #Biomass is too small for there to be any productivity
    if fOldBiomass<=minRelBmass:
        #Biomass can never be less than 0.001
        result=max([fOldBiomass,0.001])
        return(result)
    
    #Stock is in a near-virgin state, no possible increase in biomass
    if fOldBiomass>=maxRelBmass:
        return(maxRelBmass)
        
    #Actual calculations to do
    
    #Time increment
    deltaT=NewTime-OldTime
    #Maximum time-increment is one year
    if deltaT>1:
        SplitTime=[(OldTime+NewTime)/2,NewTime]
        
        #Do the increment in two steps
        SplitEst=NextRelBmass(OldBiomass=fOldBiomass,OldTime=OldTime,NewTime=SplitTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=None,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
        #Use the result corresponding to the second half of the interval
        result=SplitEst[-1]
        return(result)
    
    #Ignore time intervals less than minimum increment
    if deltaT<MinTimeInc:
        return(fOldBiomass)
    
    #Calculate the first four derivatives of relative biomass with respect to time
    coef=deriv(fOldBiomass,a,b,xmax,fmax,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
    
    #terms for the taylor polynomial
    term1= deltaT    *coef[1]
    term2=(deltaT**2)*coef[2]/2.    
    term3=(deltaT**3)*coef[3]/6.
    term4=(deltaT**4)*coef[4]/24.
    
    #Candidate estimate of change in biomass.
    deltaB=term1+term2+term3+term4
    #deltaB=max([0,deltaB])#Productivity must be positive
    
    #Candidate for new biomass
    CanEst=fOldBiomass+deltaB
    
    #Use 4th term as estimate of amount of error.
    #Split the time-interval in half if too much error occurs.  The splitting will be recursive.
    if (    (deltaB<0) or \
            abs(term4)>AllowableIntegrationError) or \
            ( (deltaB>0) and  (abs(term4)/deltaB)>AllowableIntegrationError) or \
            (deltaB>((1-fOldBiomass)/4)) or \
            (deltaB>(fOldBiomass/4)):
        SplitTime=[(OldTime+NewTime)/2,NewTime]
        
        #Do the increment in two steps
        SplitEst=NextRelBmass(OldBiomass=fOldBiomass,OldTime=OldTime,NewTime=SplitTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=None,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
        #Use the result corresponding to the second half of the interval
        result=SplitEst[-1]
        return(result)
    
    
    #Force Estimate to be in a valid range
    CanEst=ApplyBounds(CanEst,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
    
    #The candidate estimate is acceptable
    return(CanEst)
        
def deriv(CurBiomass,a,b,xmax,fmax,minRelBmass=1.e-3,maxRelBmass=1-1e-3):
    '''
    deriv(CurBiomass,a,b,xmax,fmax,minRelBmass=1.e-3,maxRelBmass=1-1e-3)
    
    Calculate first four derivatives of biomass with respect to time.
    '''   
    #results will be in a dictionary. Keys correspond to 1st, 2nd, 3rd, 4th derivative
    result={}
    
    #Derivatives are set to zero if biomass is too low or too high
    if (CurBiomass<=minRelBmass) or (CurBiomass>=maxRelBmass):
        result[1],result[2],result[3],result[4]=0., 0., 0., 0.
        return(result)
    
    #The first derivative is just the WCH productivity function
    try:
        result[1]=fmax*\
                    ((   CurBiomass /   xmax)**a) *\
                    (((1-CurBiomass)/(1-xmax))**b)    
    except RuntimeWarning:
        print('NextRelBmass 135 ', CurBiomass,fmax,a,b,xmax)
        result[1]=fmax*\
                    ((   CurBiomass /   xmax)**a) *\
                    (((1-CurBiomass)/(1-xmax))**b)
    t2=  a/CurBiomass                       -   b/(1-CurBiomass)
    t3= -a/CurBiomass/CurBiomass            -   b/(1-CurBiomass)/(1-CurBiomass)
    t4=2*a/CurBiomass/CurBiomass/CurBiomass - 2*b/(1-CurBiomass)/(1-CurBiomass)/(1-CurBiomass)
    #Other derivatives
    result[2]=result[1]*result[1]*t2
    result[3]=result[1]*result[1]*result[1]*(2*t2*t2+t3)
    result[4]=result[1]*result[1]*result[1]*result[1]*(6*t2*t2*t2+7*t2*t3+t4)
    
    return(result)            

def CurRelBmass(NewTime=None,a=None,b=None,xmax=None,fmax=None,VBiomass=None,HarvestHist=None,minRelBmass=1.e-3,maxRelBmass=1-1e-3):
    """
    CurRelBmass(NewTime=None,a=None,b=None,xmax=None,fmax=None,VBiomass=None,HarvestHist=None,minRelBmass=1.e-3,maxRelBmass=1-1e-3)
    Calculate biomass at time=NewTime given a harvest history of HarvestHist and an assumption that the stock started in a virgin state.
    
    HarvestHist is a list with the same structure as  Harvest.hdata   
        [  [decimal year, harvest amount],
           [decimal year, harvest amount],
           ...,
           [decimal year, harvest amount]]

    
    Similar to NextRelBmass except that starting biomass is fixed and multiple harvests are allowed.
    
    For single simulations, it is more efficient to use NextRelBmass
    For Bayesian models, CurRelBmass results in a much simpler relationships between nodes.  Some memory issues can be avoided.
    """
    #If there is no harvest then stock will always be in virgin state
    if HarvestHist==None:
        return(maxRelBmass)
        
    #Force time of survey to be a decimal year
    dyear=ToDecimalYear(NewTime)
        
    #Harvest history prior to survey
    useHarvestHist=[t    for t in HarvestHist if t[0]< dyear]

    #Virgin state if no harvest prior to survey        
    if not(useHarvestHist):
        return(maxRelBmass)

    #Relative Biomass after first harvest
    OldBiomass=ApplyBounds(1.-useHarvestHist[0][1]/VBiomass,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
    OldTime=useHarvestHist[0][0]
        
    #Go through the rest of the harvest records until reach survey-date
    for HarvRecord in useHarvestHist[1:]:
            
        
        #Simulate to the next harvest event
        nt=HarvRecord[0]
        NewBiomass=NextRelBmass(OldBiomass=OldBiomass,OldTime=OldTime,NewTime=nt,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=None,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
           
        #Apply harvest
        NewBiomass-=HarvRecord[1]/VBiomass
        NewBiomass=ApplyBounds(NewBiomass,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)
          
        OldBiomass=NewBiomass
        OldTime=nt

    #Continue to survey after last harvest
    NewBiomass=NextRelBmass(OldBiomass=OldBiomass,OldTime=OldTime,NewTime=NewTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBiomass,CurHarvest=None,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)  
    NewBiomass=ApplyBounds(NewBiomass,minRelBmass=minRelBmass,maxRelBmass=maxRelBmass)           
    return(NewBiomass)            
            
def ApplyBounds(x,minRelBmass=1.e-3,maxRelBmass=1-1e-3):
    result=min([x,maxRelBmass])
    result=max([result,minRelBmass])
    return(result)

if __name__ == "__main__":  
    #print('NextRelBmass 39', fOldBiomass,OldTime,NewTime,a,b,xmax,fmax,VBiomass,CurHarvest)
    
    
    a,b,xmax,fmax,VBmass_16=0.0227450862669,0.0151897315899,0.599583378857,0.385205221961,11619.5219535
    HarvestHist=CurHarvest=[[1999.0, 1378.01362006], [1999.0849315068492, 3142.94153173], [2000.0846994535518, 3476.3319236800003], [2001.0849315068492, 3420.99365454], [2002.0, 874.97968173], [2002.0849315068492, 2414.47218551], [2003.0849315068492, 2345.0725529], [2004.0846994535518, 2802.74725423], [2005.0849315068492, 3034.5329553], [2006.1616438356164, 2757.8416096]]
    
 
        
    NewTime=2000.0846994535518
    OldTime=1999.0849315068492
    OldBiomass=0.6422622219351948
    CurHarvest=0
    NewBiomass=NextRelBmass(OldBiomass=OldBiomass,OldTime=OldTime,NewTime=NewTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBmass_16,CurHarvest=CurHarvest,minRelBmass=1.e-3,maxRelBmass=1-1e-3)
    
    print()
    print(OldTime,NewTime)
    print(OldBiomass,NewBiomass)       
    
    
    
    testtime=[t for t in range(1998,2017)]    
    
    t=2000
    x=  CurRelBmass(NewTime=t,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBmass_16,HarvestHist=HarvestHist,minRelBmass=1.e-3,maxRelBmass=1-1e-3)  
    
    
    for t in testtime:
        print(t,CurRelBmass(NewTime=t,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBmass_16,HarvestHist=HarvestHist,minRelBmass=1.e-3,maxRelBmass=1-1e-3))
    
    print()
    OldBiomass=0.999
    OldTime=1998
    CurHarvest=0
    for t in   HarvestHist:  
        NewTime=t[0]
        NewBiomass=NextRelBmass(OldBiomass=OldBiomass,OldTime=OldTime,NewTime=NewTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBmass_16,CurHarvest=CurHarvest,minRelBmass=1.e-3,maxRelBmass=1-1e-3)
        print(OldTime,NewTime)
        print(OldBiomass,NewBiomass)
        print()
        OldBiomass=NewBiomass-t[1]/VBmass_16
        OldTime=NewTime
        
    OldTime=2000.0846994535518
    OldBiomass=0.6422622219351948
    NewBiomass=NextRelBmass(OldBiomass=OldBiomass,OldTime=OldTime,NewTime=NewTime,a=a,b=b,xmax=xmax,fmax=fmax,VBiomass=VBmass_16,CurHarvest=CurHarvest,minRelBmass=0,maxRelBmass=1-1e-6)
    
    print()
    print(OldTime,NewTime)
    print(OldBiomass,NewBiomass)        