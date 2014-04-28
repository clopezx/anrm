# Fits ANRM 1.0 (Irvin et. al 2013) against single-cell measurements
# of caspase reporters.

import pickle
import bayessb
import random as ra 
import numpy as np
import calibratortools as ct
import simulator_1_0 as sim
import bayes_mcmc as bmc
import matplotlib.pyplot as plt

Exp_name = ('CompII_Hyp_123_Bid_Hyp0_newtopology_1run_v0')

#----Describing Published Data-----
def ydata_fn():
    """return and array of synthesized experimental data. The array is loosely based on published experiments"""
    Apop1_td = 6.0 #six hours
    Apop2_td = 4.0 #four hours
    Necr1_td = 4.0 #four hours
    
    switchtime_CytoC = 1.0 # [hrs]
    switchtime_cPARP = 0.5 #one hour
    switchtime_MLKL = 1.0 # [hrs]
    
    Apop1_obs = ['Obs_CytoC'] #Zhang et al. Monitored CytoC (Obs_CytoC) but CytoC may not have switch behavior.
    Apop2_obs = ['Obs_cPARP']
    Necr1_obs = ['Obs_MLKL']
    
    ydata = {}
    ydata['Apop1'] = [np.array([[(Apop1_td-2*switchtime_CytoC),(Apop1_td-switchtime_CytoC), (Apop1_td-switchtime_CytoC/2), (Apop1_td-switchtime_CytoC/4), (Apop1_td-switchtime_CytoC/8), Apop1_td, (Apop1_td+switchtime_CytoC/8), (Apop1_td+switchtime_CytoC/4), (Apop1_td+switchtime_CytoC/2), (Apop1_td+switchtime_CytoC)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1],[0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Apop1_obs]
    ydata['Apop2'] = [np.array([[(Apop2_td-2*switchtime_cPARP),(Apop2_td-switchtime_cPARP), (Apop2_td-switchtime_cPARP/2), (Apop2_td-switchtime_cPARP/4), (Apop2_td-switchtime_cPARP/8), Apop2_td, (Apop2_td+switchtime_cPARP/8), (Apop2_td+switchtime_cPARP/4), (Apop2_td+switchtime_cPARP/2), (Apop2_td+switchtime_cPARP)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1], [0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Apop2_obs]
    ydata['Necr1'] = [np.array([[(Necr1_td-2*switchtime_MLKL),(Necr1_td-switchtime_MLKL), (Necr1_td-switchtime_MLKL/2), (Necr1_td-switchtime_MLKL/4), (Necr1_td-switchtime_MLKL/8), Necr1_td, (Necr1_td+switchtime_MLKL/8), (Necr1_td+switchtime_MLKL/4), (Necr1_td+switchtime_MLKL/2), (Necr1_td+switchtime_MLKL)], [0, 0, 0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95, 1], [0.025, 0.025, 0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05, 0.025]]).T, Necr1_obs]
    
    return ydata

#----Normalize--------------
ydata = ydata_fn()
ynorm = ydata.copy()
normalize = ct.normalize_array
for k in ynorm.keys():
    ynorm[k] = [normalize(ynorm[k][0], option = 1), ynorm[k][1]]


# plot data
plt.ion()
tspan = np.linspace(0,36000,1000)/3600
ii = 0
colors = ['b', 'g', 'r', 'c']

label = {}
label['Apop1'] = 'CytoC in WT Jurkat Cells with 10ng/mL TNF'
label['Apop2'] = 'cPARP in WT Jurkat Cells with 20ng/mL TNF'
label['Necr1'] = 'MLKL in -/-FADD Jurkat Cells with 30ng TNF and 20uM zVAD'

for k in ydata.keys():
    plt.errorbar(ynorm[k][0][:,0], ynorm[k][0][:,1], yerr = ynorm[k][0][:,2], fmt = '%s.' % colors[ii], label = '%s Data' % label[k])

    yinitial = pickle.load(open('%s_Initial_Values_%s.pkl' % (Exp_name, k)))
    plt.plot(tspan, yinitial, '%s--' % colors[ii], label = 'Initial for %s' % label[k])

    yfinal = pickle.load(open('%s_Final_Values_%s.pkl' % (Exp_name, k)))
    plt.plot(tspan, yfinal, '%s-' % colors[ii], label = 'Final for %s' % label[k])

    ii = ii+1

plt.xlabel('time [hrs]')
plt.title('Apoptotic and Necrotic Signals')
plt.legend(loc = 'lower left', bbox_to_anchor = (-0.1, -0.5))
