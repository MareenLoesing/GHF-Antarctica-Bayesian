import numpy as np
import toy as toy
import matplotlib.pyplot as plt
import pandas as pd
import time
plt.ion()
from Tpy import TempCalc,TempCalcExp

#Whole Antarctic continent --> Default
#Only Peninsula with heat production 
#from Burton-Johnson et al. (2017) --> Peninsula
what = 'Peninsula'
model = 'AN1'

if (what == 'Default'):
    data = pd.read_csv('Models.csv', sep='\t')
    p=4
    print('Default')
elif (what == 'Peninsula'): 
    data = pd.read_csv('Peninsula_HP.csv', sep='\t')
    A = data.HP*1e-6
    p=3
    print('Peninsula')
    
T0=0 #Surface Temperature
N = len(data)
Iterations = 1000
all_means = np.zeros((N,p))
all_std = np.zeros((N,p))
all_accept = np.zeros((N,p)) 

Temp = np.zeros((N,2))
goal_T = np.array([580.0,1315.0])
misfit = toy.IndependentGaussianMisfit(goal_T)
hyper_mean = np.zeros(N)

if (what=='Peninsula'):   
    prior = toy.RangePrior([1.0,2.0,0.0],[3.0,4.0,200e-3]) 
    proposal = toy.ComponentUpdateProposal([0.05,0.05,1e-3])
else:
    prior = toy.RangePrior([1.0,2.0,0.0,0.0],[3.0,4.0,200e-3,7.0e-6])
    proposal = toy.ComponentUpdateProposal([0.05,0.05,1e-3,1.0e-7])

hyperprior = toy.RangePrior([1.0],[10.0])
hyperProposal = toy.ComponentUpdateProposal([0.5])

for i in range(N):
    curie = data.AN1_Curie[i]*1000
    moho = data.AN1_Moho[i]*1000
    lab = data.AN1_LAB[i]*1000

    if (what == 'Peninsula'):
        Aa = A[i]
        forward = TempCalc(moho,[curie,lab],0.0,Aa) #or TempCalcExp
        x0 = np.array([2.5,3.0,50e-3])
    else: 
        forward = TempCalc(moho,[curie,lab],0.0) #or TempCalcExp
        x0 = np.array([2.5,3.0,50e-3,1.5e-6])

    
    t_0 = time.time()
    xchain,hyperChain,Lchain,accepted = toy.MCMC(x0,forward,misfit,prior,
                                                 proposal,Iterations,
                                                 np.array([5.0]),hyperProposal,
                                                 hyperprior)
    t_1 = time.time()
    Percentage = np.round(i*100/N,1)
    print ('It took %.1f seconds' % (t_1-t_0),Percentage,'% of Task complete')
    all_means[i,:] = xchain[500:,:].mean(0) 
    all_std[i,:] = xchain[500:,:].std(0) 
    hyper_mean[i] = hyperChain[500:].mean(0)
    Temp[i,:] = forward(all_means[i,:])
    print(Temp[i,0], Temp[i,1])
       
Tc = Temp[:,0]
Tl = Temp[:,1]
    
    
M = np.vstack((data.Lon,data.Lat,all_means[:,0:p].T,all_std[:,0:p].T,
               Tc,Tl,hyper_mean)).T
               
for j in range(p):
    all_accept[i,j]=np.sum(np.abs((np.diff(xchain[:,j],axis=0)))>0)*p*2.0/Iterations
print ('Probability, with which next modelparameter is taken', all_accept[i,:])
         
np.savetxt('Results_%s_%s.txt' % (model, what), M, fmt="%.3e")    



##-------------------------------Plot------------------------------------------

x=data.Lon
y=data.Lat
k1=all_means[:,0]
k2=all_means[:,1]
qD=all_means[:,2]
if (what=='Default'):
    A=all_means[:,3]

qS=1000*(qD+A*moho)
print(qS[10])


fig, axs = plt.subplots(5, sharex=True, sharey=True,figsize=(5,15))
fig.suptitle('Parameters Inversion')
plot1 = axs[0].scatter(x,y,c=k1,cmap='jet', marker=',', s=5)
axs[0].set_title('Crustal Thermal Conductivity')
cb1 = plt.colorbar(plot1,ax=axs[0])
cb1.set_label('[W/mK]', labelpad=20, rotation=270) 

plot2 = axs[1].scatter(x,y,c=k2,cmap='jet', marker=',', s=5)
axs[1].set_title('Mantle Thermal Conductivity')
cb2 = plt.colorbar(plot2,ax=axs[1])
cb2.set_label('[W/mK]', labelpad=20, rotation=270) 

plot3 = axs[2].scatter(x,y,c=qD*1000,cmap='jet', marker=',', s=5)
axs[2].set_title('Mantle Heat Flux')
cb3 = plt.colorbar(plot3,ax=axs[2])
cb3.set_label('[mW/m$^2$]', labelpad=20, rotation=270) 

plot4 = axs[3].scatter(x,y,c=A*1e6,cmap='jet', marker=',', s=5)
axs[3].set_title('Heat Production')
cb4 = plt.colorbar(plot4,ax=axs[3])
cb4.set_label('[$\mu$Wm$^{-3}$]', labelpad=20, rotation=270) 

plot5 = axs[4].scatter(x,y,c=qS,cmap='jet', marker=',', s=5, vmin=10, vmax=150)
axs[4].set_title('Surface Heat Flux')
cb5 = plt.colorbar(plot5,ax=axs[4])
cb5.set_label('[mW/m$^2$]', labelpad=20, rotation=270) 

fig.savefig('Parameter_%s_%s.png' % (model, what))



















