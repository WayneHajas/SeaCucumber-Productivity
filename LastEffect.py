'''
2017-11-09
    New distribution for the dummy data corresponding to last effect.
    This distribution, combined with sigma^2/n, will result in the same
    prior as for the other elements in the array.
2018-03-21
   Correct error in formula for loglikelihood of LastEffect-function    
    '''


from numpy import log
from pymc import  *

def LastEffect_like(value,mu,sigma,n):
    if n<2:
        return -np.inf
    LogLik=-(n-2)/(n-1)/2*(value-mu)*(value-mu)/sigma/sigma+.5*log(n-1)
    return(LogLik)
    

def LastEffect(name,value=0, mu=0,sigma=1,n=5):
	result=Stochastic(\
		logp=LastEffect_like,\
		observed=True,\
		value=value,\
		dtype=float,\
		name=name,\
		random=None,\
		parents={'mu':mu,'sigma':sigma,'n':n},\
		doc='Corrects log-likelihood for a sample-size n to what it would be for size of one')
	return(result)