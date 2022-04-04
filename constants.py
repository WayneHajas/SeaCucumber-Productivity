'''Useful constants that are give a global-scope'''
from sys import maxsize 
    
class constants():

    MinYear=1995
    MaxYear=9999 
    ConvFactor=0.45359237
    NoProdTol=1.e-6
    TimeInc=1./12.
    t0=0
    QuadArea=20
    QuadLen=5
    def __init__(self):
       return
    
global ConstantValues
ConstantValues=constants()


if __name__ == "__main__":  

   print(constants.MinYear)
   constants.MinYear=1000
   print(constants.MinYear)
