# Simulator 1.0 (Irvin et. al 2013) simulates ANRM 

import numpy as np
import calibratortools as ct



_use_pysb = False
try:
    import pysb.core
    import pysb.integrate
    _use_pysb = True
except ImportError:
    pass

"""Options for defining a bayessb.MCMC project/run.
    
    Constructor takes no options. Interface is via direct manipulation of
    attributes on instances.
    
    Attributes
    ----------
    model : pysb.Model (or similar)
    The model to estimate. If you do not wish to use a PySB model, you may
    instead provide any object with a `parameters` attribute holding a list
    of all model parameters. The parameter objects in turn must each have a
    `value` attribute containing the parameter's numerical value. If you are
    not using a PySB model you must rely on your own code to simulate the
    model in your likelihood function instead of calling `MCMC.simulate`.
    estimate_params : list of pysb.Parameter
    List of parameters to estimate, all of which must also be listed in
    `model.parameters`.
    initial_values : list of float, optional
    Starting values for parameters to estimate. If omitted, will use the
    nominal values from `model.parameters`.
    tspan : list of float
    List of time points over which to integrate the model. Ignored if not
    using a PySB model.
    start_random : bool, optional
    Whether to start from a random point in parameter space. Defaults to
    false. (NOT IMPLEMENTED)
    boundary_option : bool, optional
    Whether to enforce hard boundaries on the walk trajectory. Defaults to
    false. (NOT IMPLEMENTED)
    rtol : float or list of float, optional
    Relative tolerance for ODE solver.
    atol : float or list of float, optional
    Absolute tolerance for ODE solver.
    
    """

class Settings():
    def __init__(self):
        self.model              = None
        self.estimate_params    = None
        self.initial_values     = None
        self.tspan              = None
        self.rtol               = None
        self.atol               = None
    
    def copy(self):
        new_options = Settings()
        new_options.__dict__.update(self.__dict__)
        return new_options

class Solver(object):
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
    
    def validate(self, options):
        """Return a validated copy of options with defaults applied."""
        
        options = options.copy()
        
        if options.model is None:
            raise Exception("model not defined")
        
        if options.estimate_params is None or not len(options.estimate_params):
            raise Exception("estimate_params must contain a list of parameters")
        return options

    def run(self):
        """Initialize internal state and runs the parameter estimation."""
        self.initialize()

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
            
        # we actually work in a log-transformed phase space
        self.initial_position = np.log10(self.initial_values)
        self.position = self.initial_position
        
        # need to do this before init_solver
        self.ode_options = {};
        if self.options.rtol is not None:
            self.ode_options['rtol'] = self.options.rtol
        if self.options.atol is not None:
            self.ode_options['atol'] = self.options.atol
        
        # create solver so we can calculate the posterior
        self.init_solver()
        return self.initial_values

    def init_solver(self):
        """Initialize solver from model and tspan."""
        if _use_pysb and isinstance(self.options.model, pysb.core.Model):
            self.solver = pysb.integrate.Solver(self.options.model,
                                                self.options.tspan,
                                                **self.ode_options)

    def simulate(self, position=None, observables=False, initial_conc = None):
        """Simulate the model.
            
            Parameters
            ----------
            position : list of float, optional
            log10 of the values of the parameters being estimated. (See
            the `cur_params` method for details)
            observables : boolean, optional
            If true, return a record array containing the trajectories of the
            model's observables. If false, return a float array of all species 
            trajectories. Defaults to false.
            
            """

        if position is None:
            position = self.position

        if initial_conc is not None:
            yi = np.zeros(len(self.options.model.species))
            
            for cp, ic_param in self.options.model.initial_conditions:
                
                ii = self.options.model.parameters_initial_conditions().index(ic_param)
                si = self.options.model.get_species_index(cp)
                yi[si] = initial_conc[ii]
            self.solver.run(self.cur_params(position), yi)
        else:
            self.solver.run(self.cur_params(position))
        
        if observables:
            return self.solver.yobs
        else:
            return self.solver.y


    def cur_params(self, position=None):
        """Return a list of the values of all model parameters.
            
            For a given set of values for the parameters to be estimated, this
            method returns an array containing the actual (not log-transformed)
            values of all model parameters, not just those to be estimated, in the
            same order as specified in the model. This is helpful when simulating
            the model at a given position in parameter space.
            
            Parameters
            ----------
            position : list of float, optional
            log10 of the values of the parameters being estimated. If omitted,
            `self.position` (the most recent accepted MCMC move) will be
            used. The model's nominal values will be used for all parameters
            *not* being estimated, regardless.
        
            """
        if position is None:
            position = self.position
            # start with the original values
        values = np.array([p.value for p in self.options.model.parameters])
            # now "overlay" any rates we are estimating, by extracting them from
            # position and inverting the log transform
        values[self.estimate_idx] = np.power(10, position)
        return values

"""
    TODO
    assign all the variables that will be used in the simulator
"""