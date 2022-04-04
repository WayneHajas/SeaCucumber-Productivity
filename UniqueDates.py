# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 08:48:41 2021

@author: HajasW

Original version got lost
Function only gets used in EFAsite
I am rewriting it here

"""

from ToDecimalYear import ToDecimalYear
def UniqueDates(listDates):
    '''
     Parameters
    ----------
    listDates : TYPE
        every element is a date.

    Returns
    -------
    unique members of list. Expressed as an array of decimaly years decimal-year.

    '''
    #Make sure dates exist as decimal-years
    dyear=ToDecimalYear(listDates)
    
    #Unique values
    uniqueVal=list(set(dyear))
    
    #Sort results for good measure
    uniqueVal.sort()
    return(uniqueVal)