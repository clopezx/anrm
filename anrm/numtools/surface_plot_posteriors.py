# Fits ANRM 1.0 (Irvin et. al 2013) against single-cell measurements
# of caspase reporters.

import pickle
import bayessb.plot as bp
import random as ra 
import numpy as np
import calibratortools as ct
import simulator_1_0 as sim
import bayes_mcmc as bmc
import matplotlib.pyplot as plt

from anrm.irvin_mod_v5_tester import model

#-----------Previously Calibrated Parameters------------
initial_position = pickle.load(open('CompII_Hypthesis_123_newtopology_2run_v40_Position.pkl'))

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

#----Objective Function----------
def objective_fn(position):
    """return the value of the objective function"""
    objective = []
    for k in conditions.keys():
        ysim = solve.simulate(position, observables=True, initial_conc=conditions[k])
        PARP_MLKL_signals   = ct.extract_records(ysim, ['Obs_cPARP', 'Obs_MLKL'])
        
        if (k == 'BidKO'):
            if max(PARP_MLKL_signals[0]>0):
                td_PARP = ct.calculate_time_delay(PARP_MLKL_signals[:,0], sims.tspan)
                td_MLKL = ct.calculate_time_delay(PARP_MLKL_signals[:,1], sims.tspan)
                if td_PARP < td_MLKL:
                    objective.append(abs(td_PARP - td_MLKL))
        
        else:
            ysim_array = ct.extract_records(ysim, ydata_norm[k][1])
            ysim_norm  = ct.normalize(ysim_array, option = 1)
            ysim_tp    = ct.cubic_spline(solve.options.tspan, ysim_norm, ydata_norm[k][0][:,0]*3600)
            
            if (k == 'Necr1'):
                objective.append(np.sum((ydata_norm[k][0][:,1] - ysim_tp) ** 2 / (2 * ydata_norm[k][0][:,2])))
            
            else:
                td_PARP = ct.calculate_time_delay(PARP_MLKL_signals[:,0], sims.tspan)
                td_MLKL = ct.calculate_time_delay(PARP_MLKL_signals[:,1], sims.tspan)
                if td_MLKL < td_PARP:
                    objective.append(np.sum((ydata_norm[k][0][:,1] - ysim_tp) ** 2 / (2 * ydata_norm[k][0][:,2]))+abs(td_PARP - td_MLKL))
                else:
                    objective.append(np.sum((ydata_norm[k][0][:,1] - ysim_tp) ** 2 / (2 * ydata_norm[k][0][:,2])))
    
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

#----Normalize--------------
ydata = ydata_fn()
ydata_norm = ydata.copy()
normalize = ct.normalize_array
for k in ydata_norm.keys():
    ydata_norm[k] = [normalize(ydata_norm[k][0], option = 1), ydata_norm[k][1]]

#----Initial Protein Concetrations----
conditions = {}
init_conc = {'Apop1':{'TNFa_0': 600}, 'Apop2':{'TNFa_0': 1200}, 'Necr1':{'TNFa_0':1800, 'zVad_0':9.6e6, 'FADD_0':0}, 'BidKO':{'Bid_0': 0}}
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
opts.nsteps = 10
opts.initial_values = np.power(10, initial_position)
#opts.initial_values = solve.initial_values
opts.likelihood_fn = objective_fn
opts.prior_fn = prior
opts.step_fn = step
opts.use_hessian = False
opts.hessian_period = opts.nsteps / 10
opts.seed = 2
opts.initial_conc = conditions
opts.estimate_params = model.parameters_rules()

# values for prior calculation
prior_mean = [p.value for p in solve.options.estimate_params]
prior_ln_mean = np.log10(prior_mean)
prior_var = 6.0

mcmc = bmc.MCMC(opts)
mcmc.run()

dim0 = 0
dim1 = 1
experiment ='CompII_Hypthesis_123_newtopology_2run_v40'

# show prediction for C trajectory, which was not fit to
bp.surf(mcmc, dim0, dim1, experiment)

"""
Display the posterior of an MCMC walk on a 3-D surface.
    
Parameters
----------
mcmc : bayessb.MCMC
    MCMC object to display.
dim0, dim1 : indices of parameters to display
mask : bool/int, optional
    If True (default) the annealing phase of the walk will be discarded
    before plotting. If False, nothing will be discarded and all points will
    be plotted. If an integer, specifies the number of steps to be discarded
    from the beginning of the walk.
walk : bool, optional
    If True (default) render the walk positions. If False, do not render it.
    rejects : bool, optional
    If True (default) render each rejected position with an 'X'. If False,
    do not render them.
step : int, optional
    Render every `step`th positions along the walk. Defaults to 1 (render
    all positions). Useful to improve performance with very long walks.
square_aspect : bool, optional
    If True (default) the X and Y scales of the plot will be equal, allowing
    for direct comparison of moves in the corresponding parameter axes. If
    False the scales will auto-adjust to fit the data tightly, allowing for
    visualization of the full variance along both axes.
margin : float, optional
    Fraction of the X and Y ranges to add as padding to the surface, beyond
    the range of the points in the walk. Defaults to 0.1. Negative values
    are allowed.
    bounds0, bounds1 : array-like, optional
    Explicit ranges (min, max) for X and Y axes. Specifying either disables
    `square_aspect`.
zmin, zmax : float, optional
    Max/min height (posterior value) for the sampled surface, and the limits
    for the Z axis of the plot. Any surface points outside this range will
    not be rendered. Defaults to the actual range of posterior values from
    the walk and the sampled surface.
position_base : array-like, optional
    Vector in log10-parameter space providing values for dimensions *other*
    than dim0/dim1 when calculating the posterior surface (values at
    position dim0 and dim1 will be ignored). Defaults to the median of all
    positions in the walk.
parallelize : bool, optional
    If True (default), use the multiprocessing module to calculate the
    posterior surface in parallel using all available CPU cores. If False,
    do not parallelize.
gridsize : int, optional
    Number of points along each axis at which to sample the posterior
    surface. The total number of samples will be `gridsize`**2. Defaults to
    20. Increasing this value will produce a smoother posterior surface at
    the expense of more computational time.
    """



