"""
Apoptosis-Necrosis Reaction Network Model: 

"""

import numpy as np
import pylab as p
import calibratortools as ct
import simulator_1_0 as sim

from pysb.integrate  import odesolve
from anrm.irvin_anrm_bid_experiment_0 import model 

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,86400,1000) #24hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-3
sims.atol = 1e-6

solve = sim.Solver(sims)
solve.run()

#-----------Initial Conditions--------------------------
ic_params  = model.parameters_initial_conditions()
conditions = ct.initial_conditions(['Bak_0', 'Bax_0', 'Bid_0'], [0, 0, 0], ic_params)
ysim = solve.simulate(position = None, observables=True, initial_conc = conditions)

yout = ct.extract_records(ysim, ['Obs_cPARP', 'Obs_MLKL'])

p.ion()
p.plot(sims.tspan, yout[:,0], label = 'Cleaved Parp')
p.plot(sims.tspan, yout[:,1], label = 'MLKL')

p.xlabel('time [sec]')
p.ylabel('PARP and MLKL concentration [molecules per cell]')
p.legend()