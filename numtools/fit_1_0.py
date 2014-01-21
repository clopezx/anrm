# Fits EARM 1.0 (Albeck et. al 2008) against single-cell measurements
# of caspase reporters and MOMP timing.

import bayessb
import pysb.integrate
import numpy
import matplotlib.pyplot as plt
import os
import itertools
import scipy

from anrm.irvin_AN_crosstalk import model
from scipy.interpolate import *


def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    #return (trajectories - ymin) / (ymax - ymin)
    #normalizing the trajectories like this makes it impossible to not how scale differs
    #cell types.
    return (trajectories/ymax)
def normalize_var(trajectories, var):
    """Rescale a list of the variance of data, to match the normalized data trajectories"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return var / ymax**2

def extract_records(recarray, names):
    """Convert a record-type array and list of names into a float array"""
    return numpy.vstack([recarray[name] for name in names]).T

def likelihood(mcmc, position, initial_conc, i=None):
    """Distance between model trajectories and experimental data"""
    ysim = mcmc.simulate(position, observables=True, initial_conc=initial_conc)
    ysim_array = extract_records(ysim, obs_names)
    ysim_norm = normalize(ysim_array)
    
    ysim_objf = cubic_spline(tspan, ysim_norm, timepoints)
    print ysim_objf
    return numpy.sum((ydata_norm[:,i] - ysim_objf) ** 2 / (2 * exp_var[:,i]))

def prior(mcmc, position):
    """Distance to original parameter values"""
    return numpy.sum((position - prior_ln_mean) ** 2 / ( 2 * prior_var))

def step(mcmc):
    """Print out some statistics every 20 steps"""
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

def cubic_spline(xsim, ysim, xdata):
    """Approximates the simulated data at time = xdata using cubic spline"""
    tck = scipy.interpolate.splrep(xsim, ysim)
    yspline = scipy.interpolate.splev(xdata, tck)
    return yspline

def initial_concentrations(conditions = None):
    """provides the initial protien concentrations for the model"""
    init_conc = model.parameters_initial_conditions()
    print conditions
    init_conc_values = [p.value for p in init_conc]
    init_conc_names = [p.name for p in init_conc]
    knockout_idx = [i for i,n in enumerate(init_conc_names) if n in conditions]
    for i in knockout_idx:
        init_conc_values[i] = 0
    return init_conc_values


# data is pre-processed to filter noise and rescale the values to 0-1
data_filename1 = os.path.join(os.path.dirname(__file__), 'Zinkel_proC3_data.npy')
data_filename2 = os.path.join(os.path.dirname(__file__), 'Zinkel_proC3_var.npy')
ydata = numpy.load(data_filename1)
yvar  = numpy.load(data_filename2)
print ydata[:,1:4]
ydata_norm = normalize(ydata[:,1:4])
print ydata_norm
exp_var    = normalize_var(ydata[:,1:4], yvar[:,0:3]) #Zinkel's data is averaged with n=3.

# There are too few time points. Therefore, the solver will be inaccurate.
# Instead, I will allow the model to run 4320 time points per every 36 hours (1000 per 30000s).
timepoints = ydata[:,0]*3600 #[hrs]
t_end = timepoints.max(0)
tspan = numpy.linspace(0, t_end, 3000)

obs_names = ['Obs_pC3']

opts = bayessb.MCMCOpts()
opts.model = model
opts.tspan = tspan

opts.protien_conctrations = initial_concentrations
opts.experim_conditions = {'WT':[], 'BIDKO':['Bid_0'], 'DKO':['Bax_0', 'Bak_0'], 'TKO':['Bid_0', 'Bak_0', 'Bax_0']}

# values for prior calculation
opts.nsteps = 10000
opts.likelihood_fn = likelihood
opts.prior_fn = prior
opts.step_fn = step
opts.seed = 1
opts.atol=1e-6
opts.rtol=1e-3

scenario = 1

# A few estimation scenarios:
if scenario == 1:
    # estimate rates only (not initial conditions)
    opts.estimate_params = model.parameters_rules()
elif scenario == 2:
    # estimate first 10 rates and use hessian
    opts.estimate_params = model.parameters_rules()[0:10]
    # Warning: hessian-guidance is expensive when fitting many parameters -- the
    # time to calculate the hessian increases with the square of the number of
    # parameters to fit!
    opts.use_hessian = True
    opts.hessian_period = opts.nsteps / 3
elif scenario == 3:
    #estimate reaction parameters and initial concentration of FADD, proC8, c-Flip, TRADD and RIP1
    opts.estimate_params = model.parameters_rules()|model.parameters_initial_conditions()[2:6] |model.parameters_initial_conditions()[10:12]# |model.parameters_initial_conditions() [31:32]
elif scenario == 4:
    #estimate non-zero initial values as well as reaction perameters.
    opts.estimate_params = model.parameters_rules() | model.parameters_initial_conditions()[0:5] | model.parameters_initial_conditions()[7:11] | model.parameters_initial_conditions()[14:33]

else:
    raise RuntimeError("unknown scenario number")

# values for prior calculation
prior_mean = [p.value for p in opts.estimate_params]
prior_ln_mean = numpy.log10(prior_mean)
prior_var = 6.0

mcmc = bayessb.MCMC(opts)

mcmc.run()

# print some information about the maximum-likelihood estimate parameter set
print
print '%-10s %-12s %-12s %s' % ('parameter', 'actual', 'fitted', 'log10(fit/actual)')
fitted_values = mcmc.cur_params()[mcmc.estimate_idx]
for param, new_value in zip(opts.estimate_params, fitted_values):
    change = numpy.log10(new_value / param.value)
    values = (param.name, param.value, new_value, change)
    print '%-10s %-12.2g %-12.2g %-+6.2f' % values

#generates final and initial output values from the model and plot them.
initial_conc_WT = initial_concentrations(opts.experim_conditions['WT'])
initial_conc_BIDKO = initial_concentrations(opts.experim_conditions['BIDKO'])
initial_conc_DKO = initial_concentrations(opts.experim_conditions['DKO'])
initial_params = [p.value for p in opts.estimate_params]
yfinal_array_WT = normalize(extract_records(mcmc.simulate(observables=True, initial_conc=initial_conc_WT), obs_names))
yfinal_array_BIDKO = normalize(extract_records(mcmc.simulate(observables=True, initial_conc=initial_conc_BIDKO), obs_names))
yfinal_array_DKO = normalize(extract_records(mcmc.simulate(observables=True, initial_conc=initial_conc_DKO), obs_names))
ysim_initial_WT = mcmc.simulate(numpy.log10(initial_params), observables=True, initial_conc=initial_conc_WT)
ysim_initial_BIDKO = mcmc.simulate(numpy.log10(initial_params), observables=True, initial_conc=initial_conc_BIDKO)
ysim_initial_DKO = mcmc.simulate(numpy.log10(initial_params), observables=True, initial_conc=initial_conc_DKO)

yinitial_array_WT = normalize(extract_records(ysim_initial_WT, obs_names))
yinitial_array_BIDKO = normalize(extract_records(ysim_initial_BIDKO, obs_names))
yinitial_array_DKO = normalize(extract_records(ysim_initial_DKO, obs_names))


plt.ion()
tspan_hrs = tspan/3600
time_hrs  = timepoints/3600
data_lines_WT = plt.errorbar(time_hrs, ydata_norm[:,0], exp_var[:,0], label = 'WT data')
data_lines_BIDKO = plt.errorbar(time_hrs, ydata_norm[:,1], exp_var[:,1], label = 'BIDKO data')
data_lines_DKO = plt.errorbar(time_hrs, ydata_norm[:,2], exp_var[:,2], label = 'DKO data')
data_lines_DKO = plt.errorbar(time_hrs, ydata_norm[:,2], exp_var[:,2], label = 'TKO data')


initial_lines_WT = plt.plot(tspan_hrs, yinitial_array_WT, 'b--', label = 'initial WT')
initial_lines_BIDKO = plt.plot(tspan_hrs, yinitial_array_BIDKO, 'g--', label = 'initial BIDKO')
initial_lines_DKO = plt.plot(tspan_hrs, yinitial_array_DKO, 'r--', label = 'initial DKO')
initial_lines_DKO = plt.plot(tspan_hrs, yinitial_array_DKO, 'c--', label = 'initial TKO')


final_lines_WT = plt.plot(tspan_hrs, yfinal_array_WT, 'b', label = 'final WT')
final_lines_BIDKO = plt.plot(tspan_hrs, yfinal_array_BIDKO, 'g', label = 'final BIDKO')
final_lines_DKO = plt.plot(tspan_hrs, yfinal_array_DKO, 'r', label = 'final DKO')
final_lines_DKO = plt.plot(tspan_hrs, yfinal_array_DKO, 'c', label = 'final TKO')

plt.xlabel('time [hrs]')
plt.title('Procaspase 3 cleavage vs. time')
plt.legend(loc = 'lower left')

# plot data and simulated cleaved PARP trajectories before and after the fit
colors = ('r', 'g', 'b')
patterns = ('k--', 'k', 'k:x')





