# Written by Wolfgang Szwillus

import numpy as np
from scipy.spatial.distance import cdist
_LOG_2PI = np.log(2 * np.pi)


class ComponentUpdateProposal:
    """
    Objects of this class can be called to generate a proposal for MCMC. 
    A component is selected randomly and altered randomly 
    using a normal distribution. 
    The scale of the change is determined by class parameter scale
    """  
    def __init__(self,scales):
        self.scales = scales
        self.n = len(scales)
    def __call__(self):
        i = np.random.randint(0,self.n)
        return i,np.random.randn() * self.scales[i]

        
class IndependentGaussianMisfit:
    """
    Callable, returns logpdf
    Input:
    data: Corresponds to mean of Gaussian distribution
    hyper: Corresponds to std of Gaussian distribution
    """
    def __init__(self,data):
        self.data = data
        self.n = len(data)
    def __call__(self,y,hyper):
        delta = y-self.data
        return -self.n*np.log(hyper[0])-0.5*np.sum(delta**2/hyper[0]**2)
        #probaility that model is within sigma of misfitfunction (Gauss)
    
class RangePrior:
    def __init__(self,lb,ub):
        self.lb = lb
        self.ub = ub
    def __call__(self,x):
        if np.all((x>=self.lb) & (x<=self.ub)):
            return 1.0
        else:
            return -np.inf
        
        
def MCMC(x0,forward,misfit,prior,proposal,nruns,
            hyper0=None,hyperProposal=None,hyperPrior=None,
            chainSave=1):
    x = x0.copy()
    y = forward(x)
    if not hyper0 is None:
        hyper = hyper0.copy()
    else:
        hyper = None
    if not hyper0 is None:
        oldL = misfit(y,hyper) + prior(x) + hyperPrior(hyper)
    else:
        oldL = misfit(y) + prior(x)
    xchain = np.zeros((nruns//chainSave,x.shape[0]))
    Lchain = np.zeros((nruns//chainSave))
    xchain[0,:] = x0
    if not hyper is None:
        hyperchain = np.zeros((nruns//chainSave,hyper.shape[0]))
        hyperchain[0,:] = hyper  
    accepted = 0
    Lchain[0] = oldL
    
    for i in range(nruns-1):
        if hyper is None or np.random.rand() > 0.5:
            changedParam,dx = proposal()
            xcan = x.copy()
            xcan[changedParam] = xcan[changedParam] + dx
            ycan = forward(xcan)
            if not hyper is None:
                hypercan = hyper.copy()
        else:
            hyperChanged,dh = hyperProposal()
            hypercan = hyper.copy()
            hypercan[hyperChanged] = hypercan[hyperChanged] + dh
            xcan = x.copy()
            ycan = y.copy()
        if not hyper is None:
            newL = misfit(ycan,hypercan) + prior(xcan) + hyperPrior(hypercan)
        else:
            newL = misfit(ycan) + prior(xcan)
        acceptance = np.exp(newL - oldL)
        if np.random.rand() < acceptance:
            accepted = accepted + 1
            x = xcan.copy()
            if not hyper is None:
                hyper = hypercan.copy()
            y = ycan.copy()
            oldL = newL
        
        if (i+1)%chainSave==0:
            xchain[(i+1)//chainSave,:] = x.copy()
            if not hyper is None:
                hyperchain[(i+1)//chainSave,:] = hyper.copy()
            Lchain[(i+1)//chainSave] = oldL
    if not hyper is None:
        return xchain,hyperchain,Lchain,accepted
    else:
        return xchain,Lchain,accepted

        
