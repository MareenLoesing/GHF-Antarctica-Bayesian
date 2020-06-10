import sys
import numpy as np
sys.path.append("Y:")
import toy as toy
import matplotlib.pyplot as plt
plt.ion()
import seaborn as sns
import pandas as pd

class TempCalc:
    def __init__(self,moho,isotherm_depths,surface_temp,Heat_prod=None):
        self.isotherm_depths = isotherm_depths
        self.moho = moho
        self.T0 = surface_temp
        if not Heat_prod is None: self.A = Heat_prod
        else: self.A = 0        
    def __call__(self,args):
        args = np.asarray(args)
        '''
        k1 = crustal therm. condtuctivity
        k2 = mantle therm. conduct.
        A = crustal heat production
        qD = heatflux from mantle to crust
        '''
        if (self.A == 0):
            k1,k2,qD,A = args.T           
        else:
            k1,k2,qD = args.T
            A = self.A
        T = np.zeros(k1.shape+(len(self.isotherm_depths),))
        moho = self.moho 
        for i,z in enumerate(self.isotherm_depths):
            if z<=moho:
                T[...,i]=self.T0+(qD+A*moho)/k1*z-A/(2.0*k1)*z**2
            else:
                T[...,i]=self.T0+moho*(qD/k1+0.5*A*moho/k1-qD/k2)+qD/k2*z
        return T



class TempCalcExp:
    def __init__(self,moho,isotherm_depths,surface_temp,Heat_prod=None):
        self.isotherm_depths = isotherm_depths
        self.moho = moho
        self.T0 = surface_temp
        if not Heat_prod is None: self.A = Heat_prod
        else: self.A = 0        
    def __call__(self,args):
        args = np.asarray(args)
        if (self.A == 0):
            k1,k2,qD,A = args.T           
        else:
            k1,k2,qD = args.T
            A = self.A
        H0 = A #surface heat production
        hr = 8000.0 #scale depth
        T = np.zeros(k1.shape+(len(self.isotherm_depths),))
        moho = self.moho 
        for i,z in enumerate(self.isotherm_depths):
            if z<=moho:
                exp = np.exp(-z/hr)
                T[...,i] = self.T0 + qD*z/k1 + H0*hr**2/k1*(1-exp)
            else:
                exp1 = np.exp(-moho/hr)                
                Tm = self.T0 + qD*moho/k1 + H0*hr**2/k1*(1-exp1)
                T[...,i] = Tm + (z-moho)*qD/k2
        return T



#Test, runs when this script is executed individually
if __name__ == '__main__':
    print('Main')
    forward = TempCalc(20e3,[10e3,30e3],0.0)
      
    fake_T = forward((2.5,3.0,110e-3,2.5e-6))
    misfit = toy.IndependentGaussianMisfit(fake_T)
    
    prior = toy.RangePrior([2.0,2.5,0.0,0.0],[3.0,3.5,200e-3,3.0e-6])
    proposal = toy.ComponentUpdateProposal([0.04,0.15,2e-3,2.0e-7])

    hyperPrior = toy.RangePrior([1.0],[100.0])
    hyperProposal = toy.ComponentUpdateProposal([0.5])
    
    x0 = np.array([2.75,3.25,100e-3,1.5e-6])

    #Inversion with Markov-Chain-Monte-Carlo 
    xchain,hyperChain,Lchain,accepted = toy.MCMC(x0,forward,misfit,prior,
                                                 proposal,100000,
                                                 np.array([5.0]),hyperProposal,
                                                 hyperPrior)       
    #Probabilities with which next model is taken
    for i in range(4):
        print (np.sum(np.abs((np.diff(xchain[:,i],axis=0)))>0)*8.0/100000.0)
    print(accepted)
    
    #plots correlations and histograms
    sns.pairplot(data=pd.DataFrame(xchain[50000:]))
        
    
