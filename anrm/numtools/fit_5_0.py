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

from anrm.irvin_mod_v4_tester import model

#----Experimental Data----
"""
    ydata: dict keys = name of the experimental conditions. The items in the dict are 1x2 lists
    the first item in the list is an array of data.
    Array[:,0,:] = timepoints
    array[:,1,:] = values
    array[:,2,:] = variance
    the 3rd dimension is the observables.
    the second item in the list is the observables.
    
    init_conc: dict of experimental conditions(keys) and initial mononmer concentrations (items)
    objective_fn:
    prior:
    step:
"""
#-----------Previously Calibrated Parameters------------
initial_position = pickle.load(open('CompII_Hypthesis_123_Position.pkl'))

#----User Defined Functions-----
def ydata_fn():
    """return and array of synthesized experimental data. The array is loosely based on published experiments"""
    Apop1_td = 6.0 #six hours
    Apop2_td = 4.0 #four hours
    Necr1_td = 4.0 #four hours

    switchtime_CytoC = 0.5 # [hrs]
    switchtime_cPARP = 0.5 #one hour
    switchtime_MLKL = 1.0 # [hrs]

    Apop1_obs = ['Obs_CytoC'] #Zhang et al. Monitored CytoC (Obs_CytoC) but CytoC does not have switch behavior.
    Apop2_obs = ['Obs_cPARP']
    Necr1_obs = ['Obs_MLKL']
    
    ydata = {}
    ydata['Apop1'] = [np.array([[(Apop1_td-switchtime_CytoC/2), (Apop1_td-switchtime_CytoC/4), (Apop1_td-switchtime_CytoC/8), Apop1_td, (Apop1_td+switchtime_CytoC/8), (Apop1_td+switchtime_CytoC/4), (Apop1_td+switchtime_CytoC/2)], [0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95],[0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05]]).T, Apop1_obs]
    ydata['Apop2'] = [np.array([[(Apop2_td-switchtime_cPARP/2), (Apop2_td-switchtime_cPARP/4), (Apop2_td-switchtime_cPARP/8), Apop2_td, (Apop2_td+switchtime_cPARP/8), (Apop2_td+switchtime_cPARP/4), (Apop2_td+switchtime_cPARP/2)], [0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95], [0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05]]).T, Apop2_obs]
    ydata['Necr1'] = [np.array([[(Necr1_td-switchtime_MLKL/2), (Necr1_td-switchtime_MLKL/4), (Necr1_td-switchtime_MLKL/8), Necr1_td, (Necr1_td+switchtime_MLKL/8), (Necr1_td+switchtime_MLKL/4), (Necr1_td+switchtime_MLKL/2)], [0.05, 0.205, 0.340, 0.5, 0.659, 0.794, 0.95], [0.025, 0.05, 0.075, 0.1, 0.085, 0.065, 0.05]]).T, Necr1_obs]
    
    return ydata

def objective_fn(position):
    """return the value of the objective function"""
    objective = []
    for k in conditions.keys():
        ysim = solve.simulate(position, observables=True, initial_conc=conditions[k])
        ysim_array = ct.extract_records(ysim, ynorm[k][1])
        ysim_norm  = ct.normalize(ysim_array, option = 1)
        ysim_tp    = ct.cubic_spline(solve.options.tspan, ysim_norm, ynorm[k][0][:,0]*3600)
        
        if (k == 'Necr1'):
            objective.append(np.sum((ynorm[k][0][:,1] - ysim_tp) ** 2 / (2 * ynorm[k][0][:,2])))
        
        else:
            PARP_MLKL_signals   = ct.extract_records(ysim, ['Obs_cPARP', 'Obs_MLKL'])
        
            td_PARP = calculate_time_delay(PARP_MLKL_signals[:,0])
            td_MLKL = calculate_time_delay(PARP_MLKL_signals[:,1])
            
            if td_MLKL < td_PARP:
                objective.append(np.sum((ynorm[k][0][:,1] - ysim_tp) ** 2 / (2 * ynorm[k][0][:,2]))+abs(td_PARP - td_MLKL))
            else:
                objective.append(np.sum((ynorm[k][0][:,1] - ysim_tp) ** 2 / (2 * ynorm[k][0][:,2])))

    return np.sum(objective)
    

def calculate_time_delay(signal):
    if np.isnan(np.sum(signal)):
        return None
    else:
        norm_signal = ct.normalize(signal, option = 0)
        norm_signal = norm_signal.tolist()
        idx         = norm_signal.index(min(norm_signal, key = lambda x: abs(x-0.5)))
        return ct.cubic_spline(norm_signal[idx-3:idx+3], solve.options.tspan[idx-3:idx+3], [0.5], degree = 1)
    
def prior(mcmc, position):
    """Distance to original parameter values"""
    
    return np.sum((position - prior_ln_mean) ** 2 / ( 2 * prior_var))


def step(mcmc):
    """Print out some statistics every 20 steps"""
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

#----Experiment Name--------
Exp_name = ('CompII_Hypthesis_123_apoponly')

#----Data and conditions----
ydata = ydata_fn()
#init_conc = {'Apop1':{'TNFa_0': 600}}
init_conc = {'Apop1':{'TNFa_0': 600}, 'Apop2':{'TNFa_0': 1200}, 'Necr1':{'TNFa_0':1800, 'zVad_0':9.6e6, 'FADD_0':0}} #600 = 10ng/ml TNFa, 9.6e6 = 20uM
#init_conc = {'Apop2':{'TNFa_0': 1200}, 'Necr1':{'TNFa_0':1800, 'zVad_0':9.6e6, 'FADD_0':0}}
init_conc = {'Apop1':{'TNFa_0': 600}, 'Apop2':{'TNFa_0': 1200}}

#----Normalize--------------
ynorm = ydata.copy()
normalize = ct.normalize_array
for k in ynorm.keys():
    ynorm[k] = [normalize(ynorm[k][0], option = 1), ynorm[k][1]]

#----Initial Protein Concetrations----
conditions = {}
ic_params  = model.parameters_initial_conditions()
for k in init_conc.keys():
    conditions[k] = ct.initial_conditions(init_conc[k].keys(), init_conc[k].values(), ic_params)

#----Simulator Settings----
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,36000,1000) #10hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-5
sims.atol = 1e-5

solve = sim.Solver(sims)
solve.run()



#----Bayesian and MCMC Options----
opts = bmc.MCMCOpts()
opts.nsteps = 2000
opts.likelihood_fn = objective_fn
opts.prior_fn = prior
opts.step_fn = step
opts.seed = ra.randint(0,1000)
#opts.initial_values = np.power(10, initial_position)
opts.initial_values = solve.initial_values
opts.initial_conc = conditions

# values for prior calculation
prior_mean = [p.value for p in solve.options.estimate_params]
prior_ln_mean = np.log10(prior_mean)
prior_var = 6.0

mcmc = bmc.MCMC(opts)
mcmc.run()

# print some information about the maximum-likelihood estimate parameter set
print
print '%-10s %-12s %-12s %s' % ('parameter', 'actual', 'fitted', 'log10(fit/actual)')
fitted_values = solve.cur_params(mcmc.position)[solve.estimate_idx]
for param, new_value in zip(sims.estimate_params, fitted_values):
    change = np.log10(new_value / param.value)
    values = (param.name, param.value, new_value, change)
    print '%-10s %-12.2g %-12.2g %-+6.2f' % values

# save data
initial_params = [p.value for p in sims.estimate_params]

for k in conditions.keys():
    yinitial = ct.normalize(ct.extract_records(solve.simulate(np.log10(initial_params), observables = True, initial_conc = conditions[k]), ynorm[k][1]), option = 1)
    pickle.dump(yinitial, open('%s_Initial_Values_%s.pkl' % (Exp_name, k), 'wb'))
    
    yfinal = ct.normalize(ct.extract_records(solve.simulate(mcmc.position, observables=True, initial_conc=conditions[k]),ynorm[k][1]), option = 1)
    pickle.dump(yfinal, open('%s_Final_Values_%s.pkl' % (Exp_name, k), 'wb'))

pickle.dump(mcmc.position, open('%s_Position.pkl' % Exp_name, 'wb'))


# plot data
plt.ion()
tspan = sims.tspan/3600
initial_params = [p.value for p in sims.estimate_params]
ii = 0
colors = ['b', 'g', 'r', 'c']

"""
for k in conditions.keys():
    plt.errorbar(ynorm[k][0][:,0], ynorm[k][0][:,1], yerr = ynorm[k][0][:,2], fmt = '%s.' % colors[ii], label = '%s data' % k)

    yinitial = ct.normalize(ct.extract_records(solve.simulate(np.log10(initial_params), observables = True, initial_conc = conditions[k]), ynorm[k][1]), option = 1)
    plt.plot(tspan, yinitial, '%s--' % colors[ii], label = 'initial %s' % k)

    yfinal = ct.normalize(ct.extract_records(solve.simulate(mcmc.position, observables=True, initial_conc=conditions[k]),ynorm[k][1]), option = 1)
    plt.plot(tspan, yfinal, '%s-' % colors[ii], label = 'final %s' % k)

    ii = ii+1

plt.xlabel('time [hrs]')
plt.title('CompII Hypotheses 2; Apoptotic and Necrotic Signals')
plt.legend(loc = 'lower left', bbox_to_anchor = (1.0, -0.02))
"""

"""
TODO

In sims., check that tspan.max > timepoints.max

"""