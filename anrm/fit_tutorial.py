# The purpose of this file is to understand what is going on in fit_1_0.py and the bayessb files it calls.
# It is certainly a chance to learn object oriented programming.

import os
import numpy as np 
from irvin_mod_1_0 import model

data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data.npy')
# this creates a filename. os.path.dirname(__file__), returns the path from python's
# current directory to this file (i.e. fit_tutorial.py). For example if fit_tutorial
# were imported from Research/anrm as anrm.fit_tutorial then
# os.path.dirname(__file__) would be anrm/

# os.path.join(connects os.path.dirname(__file__) to the string variable 'exper...
# and converts it to a path/file. e.g. anrm/experimental_data.npy

# for this to work, experimental_data.npy and fit_tutorial need to be in the same
# folder.

ydata_norm = np.load(data_filename) #the .npy extension is a numpy array.
exp_var = 0.2  # take variance of experimental measurements to be constant

t_end = 5.5 * 3600  # 5.5 hours, in seconds
tspan = np.linspace(0, t_end, len(ydata_norm)) #creates a list of numbers between
# zero and t_end. There are len(ydata_norm) numbers equally spaced.
# NOTE: I'll probably want to upload a time vector so that they can be spaced by how
# they were acquired.

### From bayessb.__init__.py and is instantiated in fit_1_0 via opts = bayessb.MCMCOpts()
class MCMCOpts(object):
# object == a class here. this link gives more information
# http://docs.python.org/release/2.2.3/whatsnew/sect-rellinks.html
    # whenever MCMCOpts is instantiated, __init__ performs the functionsimmediately.
    # self is a place holder. When an object is created, it's name (e.g. "opts")
    # replaces self...
    # so.. self.model becomes opts.model.
    def __init__(self):
        self.model              = None
        self.estimate_params    = None
        self.initial_values     = None
        self.tspan              = None
        self.step_fn            = None
        self.likelihood_fn      = None
        self.prior_fn           = None
        self.nsteps             = None
        self.use_hessian        = False
        self.start_random       = False
        self.boundary_option    = False
        self.rtol               = None
        self.atol               = None
        self.norm_step_size     = 0.75
        self.hessian_period     = 25000
        self.hessian_scale      = 0.085
        self.sigma_adj_interval = None
        self.anneal_length      = None
        self.T_init             = 10
        self.accept_rate_target = 0.3
        self.sigma_max          = 1
        self.sigma_min          = 0.25
        self.sigma_step         = 0.125
        self.thermo_temp        = 1
        self.seed               = None

    def copy(self):
        # this copies the object. This does not happen immediately. You have to type
        # x = opt.copy(), and x will be a copy of opt.
        new_options = MCMCOpts()
        new_options.__dict__.update(self.__dict__)
        return new_options

opts = MCMCOpts() #creates the object, opts. this object has the list of variables(?)
# assigned in def __init__(self):
opts.model = model #reassigns opts.model to be the model imported from... instead of None
opts.tspan = tspan

opts.estimate_params = model.parameters_rules() # makes opts.esti... a dictionary with
# parameter names and Parameter('...', value=..  ) commands.
prior_mean = [p.value for p in opts.estimate_params] #creates a list of numbers param.. values.
prior_var = 6.0

opts.nsteps = 1000 # number of MCMC iterations to perform.

def likelihood(mcmc, position):
    """Distance between model trajectories and experimental data"""
    ysim = mcmc.simulate(position, observables=True)
    ysim_array = extract_records(ysim, obs_names)
    ysim_norm = normalize(ysim_array)
    return np.sum((ydata_norm - ysim_norm) ** 2 / (2 * exp_var ** 2))

def prior(mcmc, position):
    """Distance to original parameter values"""
    return np.sum((position - prior_mean) ** 2 / ( 2 * prior_var))

def step(mcmc):
    """Print out some statistics every 20 steps"""
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

opts.prior_fn = prior # this takes the function, prior(mcmc, position) and assigns it as a
# method in the object, opts
opts.likelihood_fn = likelihood #this takes the function, likelihood(mcmc, position) and puts it
# into the opts object.
opts.step_fn = step # same as prior and likelihood.
opts.seed = 1 # seed for random generator
opts.atol=1e-6 # absolulte tolerance
opts.rtol=1e-3 # relative tolerance


class MCMC(object):
    def __init__(self, options):
        self.options = self.validate(options)
    
    def __getstate__(self):
        # clear solver since it causes problems with pickling
        state = self.__dict__.copy()
        del state['solver']
        return state
    
    def __setstate__(self, state):
        # re-init the solver which we didn't pickle
        self.__dict__.update(state)
        self.init_solver()
    
    def run(self):
        """Initialize internal state and runs the parameter estimation."""
        self.initialize()
        self.estimate()
    
    def validate(self, options):
        """Return a validated copy of options with defaults applied."""
        # FIXME should this live in MCMCOpts?
        
        options = options.copy()
        
        if options.model is None:
            raise Exception("model not defined")
        
        if options.estimate_params is None or not len(options.estimate_params):
            raise Exception("estimate_params must contain a list of parameters")
        
        # clamp hessian_period to actual number of steps
        if options.use_hessian:
            options.hessian_period = min(options.hessian_period, options.nsteps)
        else:
            options.hessian_period = np.inf
        
        if options.anneal_length is None:
            # default for anneal_length if unspecified
            if options.use_hessian:
                # if using hessian, anneal until we start using it
                options.anneal_length = options.hessian_period
            else:
                # otherwise, anneal for 10% of the run
                options.anneal_length = np.floor(options.nsteps * 0.10)
        else:
            # clamp it to actual number of steps
            options.anneal_length = min(options.anneal_length, options.nsteps)
        
        # default for sigma_adj_interval if unspecified
        if options.sigma_adj_interval is None:
            # default to 10 adjustments throughout the annealing phase
            options.sigma_adj_interval = max(int(options.anneal_length / 10), 1)
        
        return options
    
    def initialize(self):
        """Initialize internal state from the option set."""
        
        # create list of starting values from initial parameter values given by
        # user. vector only contains values which are to be estimated!
        self.num_estimate = len(self.options.estimate_params)
        if self.options.initial_values is not None:
            self.initial_values = self.options.initial_values
        else:
            # if no explicit values given, take values from model
            self.initial_values = [p.value
                                   for p in self.options.estimate_params]
        # indices of parameters to be estimated
        self.estimate_idx = [i for i, p
                             in enumerate(self.options.model.parameters)
                             if p in self.options.estimate_params]
        print self.initial_values


mcmc = MCMC(opts) #where all the fun begins... MCMC is a subclass of the object opts.
# It inherits the attributes and etc from the object opts.
