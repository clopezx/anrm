"""
Apoptosis-Necrosis Reaction Network Model: 

"""

import numpy as np
import pylab as p
import calibratortools as ct
import simulator_1_0 as sim

from pysb.integrate  import odesolve
from anrm.irvin_anrm_bid_experiment_3 import model

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,10000,100) #24hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-4
sims.atol = 1e-8

solve = sim.Solver(sims)
solve.run()

#-----------Initial Conditions--------------------------
ic_params  = model.parameters_initial_conditions()
conditions = ct.initial_conditions([], [], ic_params)
ysim = solve.simulate(position = None, observables=True, initial_conc = conditions)

yout = ct.extract_records(ysim, ['Obs_cPARP', 'Obs_MLKL'])

PARP_MLKL_signals   = yout
td_PARP = ct.calculate_time_delay(PARP_MLKL_signals[:,0], sims.tspan)
td_MLKL = ct.calculate_time_delay(PARP_MLKL_signals[:,1], sims.tspan)
print td_PARP, td_MLKL

p.ion()
p.plot(sims.tspan, yout[:,0], label = 'Cleaved Parp')
p.plot(sims.tspan, yout[:,1], label = 'MLKL')

p.xlabel('time [sec]')
p.ylabel('PARP and MLKL concentration [molecules per cell]')
p.legend()