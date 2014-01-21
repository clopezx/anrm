"""
Apoptosis-Necrosis Reaction Network Model: 

"""

import pysb
import numpy as np
import pylab as p
import pickle

from anrm import merge
from pysb.integrate  import odesolve

from anrm.irvin_anrm_experiment_10 import model

#Edit parameters 
new_params = pickle.load(open('TNFa_H2_Calibrated_Params_exp10.pkl'))
param_edits = merge.Edit_Parameters(model.parameters, new_params)
merged_parameters = param_edits.merged_parameters

pysb.core.Model('m')

m.monomers = model.monomers
m.compartments = model.compartments
m.parameters = merged_parameters
m.rules = model.rules
m.observables = model.observables
m.initial_conditions = model.initial_conditions

t = np.linspace(0,20000,100)
yout = odesolve(m, t)

p.ion()
p.plot(t, yout['Obs_cPARP'], label = 'Cleaved Parp')
p.plot(t, yout['Obs_MLKL'], label = 'MLKL')

p.xlabel('time [sec]')
p.ylabel('PARP concentration [molecules per cell]')
p.legend()