# Fits ANRM 1.0 (Irvin et. al 2013) against single-cell measurements
# of caspase reporters.

import pickle
import bayessb
import numpy as np
import calibratortools as ct
import simulator_1_0 as sim
import bayes_mcmc as bmc
import matplotlib.pyplot as plt

from anrm.irvin_anrm_experiment_16 import model

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

ydata = pickle.load(open('Zinkel_Data.pkl'))
init_conc = pickle.load(open('Zinkel_Conditions.pkl'))

#----User Defined Functions-----
def objective_fn(position):
    """return the value of the objective function"""
    ysim_norm = []
    for k in conditions.keys():
        ysim = solve.simulate(position, observables=True, initial_conc=conditions[k])
        ysim_array   = ct.extract_records(ysim, ynorm[k][1])
        ysim_norm.append(ct.normalize(ysim_array, option = 1))

    objective = []
    timepoints = [0, 21600, 43200, 86400]
    keys = conditions.keys()
    for i in range(len(ysim_norm)):
        ysim_tp = ct.cubic_spline(solve.options.tspan, ysim_norm[i], timepoints) #ysim(t = expe. timepoints)

        # Modifiying the objective function to penalize, when an apoptotic cell behaves necrotically, and vice versa. I consider a cell apoptotic when the time-delay (i.e. the time at which the signal is at 50% maximum) of the apoptotic signal is less than that of the necrotic signal. This aproach works so long as the gain of the signals is not as senstive to conditions at the time-delay.
    
        PARP_MLKL_signals   = ct.extract_records(ysim, ['Obs_cPARP', 'Obs_MLKL'])
        td_PARP = calculate_time_delay(PARP_MLKL_signals[:,0])
        td_MLKL = calculate_time_delay(PARP_MLKL_signals[:,1])

        if ((keys[i] == 'TKO') & (td_PARP < td_MLKL)):# | ((keys[i] == 'WT') & (td_PARP > td_MLKL)):
            objective.append(np.sum((ynorm[keys[i]][0][:,1] - ysim_tp) ** 2 / (2 * ynorm[keys[i]][0][:,2]))+abs(td_PARP - td_MLKL))
        else:
            objective.append(np.sum((ynorm[keys[i]][0][:,1] - ysim_tp) ** 2 / (2 * ynorm[keys[i]][0][:,2])))

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

#----Normalize Data----
ynorm = ydata.copy()
normalize = ct.normalize_array
for k in ynorm.keys():
    ynorm[k] = [normalize(ynorm[k][0], option = 1), ynorm[k][1]]
del ynorm['BidKO']
del ynorm['DKO']

#----Initial Protein Concetrations----
conditions = {}
ic_params  = model.parameters_initial_conditions()
for k in init_conc.keys():
    conditions[k] = ct.initial_conditions(init_conc[k].keys(), init_conc[k].values(), ic_params)
del conditions['BidKO']
del conditions['DKO']

#----Simulator Settings----
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,86400,1000) #24hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-3
sims.atol = 1e-6

solve = sim.Solver(sims)
solve.run()



#----Bayesian and MCMC Options----
opts = bmc.MCMCOpts()
opts.nsteps = 10000
opts.likelihood_fn = objective_fn
opts.prior_fn = prior
opts.step_fn = step
opts.seed = 1
opts.initial_values = solve.initial_values
opts.initial_conc = conditions

# values for prior calculation
prior_mean = [p.value for p in solve.options.estimate_params]
prior_ln_mean = np.log10(prior_mean)
prior_var = 7.0

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

#plot data
plt.ion()
tspan = sims.tspan/3600
initial_params = [p.value for p in sims.estimate_params]
ii = 0
colors = ['b', 'g', 'r', 'c']
for k in ['WT', 'TKO']:
    plt.errorbar(ynorm[k][0][:,0], ynorm[k][0][:,1], yerr = ynorm[k][0][:,2], fmt = '%s.' % colors[ii], label = '%s data' % k)

    yinitial = ct.normalize(ct.extract_records(solve.simulate(np.log10(initial_params), observables = True, initial_conc = conditions[k]), ynorm[k][1]), option = 1)
    plt.plot(tspan, yinitial, '%s--' % colors[ii], label = 'initial %s' % k)

    yfinal = ct.normalize(ct.extract_records(solve.simulate(mcmc.position, observables=True, initial_conc=conditions[k]),ynorm[k][1]), option = 1)
    plt.plot(tspan, yfinal, '%s-' % colors[ii], label = 'final %s' % k)

    ii = ii+1

plt.xlabel('time [hrs]')
plt.title('Procaspase 3 cleavage vs. time')
plt.legend(loc = 'lower left', bbox_to_anchor = (1.0, -0.02))


"""
TODO

In sims., check that tspan.max > timepoints.max

"""