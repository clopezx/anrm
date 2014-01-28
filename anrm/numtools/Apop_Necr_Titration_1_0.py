# Runs ANRM 1.0 (Irvin et. al 2013) under a range of pro-apoptotic and pro-necrotic caspase 8
# concentrations.

import numpy as np
import pylab as pl
import calibratortools as ct
import simulator_1_0 as sim

# ----------Model and Initial Conditions----------------
from anrm.irvin_anrm_experiment_16 import model

range_proC8 = np.linspace(0,50000,51) #starting at zero causes NaNs when you normalize the data.
range_cFlip = np.linspace(0,200000,201)
range_Bid   = np.linspace(0,16000,51)
range_XIAP  = np.linspace(0,35000, 101)

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,86400,1000) #24hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-3
sims.atol = 1e-6

solve = sim.Solver(sims)
solve.run()

delta_td = []
apopt_td = []
necro_td = []

for i in range_XIAP:
    #-----------Initial Conditions--------------------------
    ic_params  = model.parameters_initial_conditions()
    conditions = ct.initial_conditions(['XIAP_0'], [i], ic_params)
    ysim = solve.simulate(position = None, observables=True, initial_conc = conditions)

    #-----------Calculate Time Delays-----------------------
    PARP_MLKL_signals   = ct.extract_records(ysim, ['Obs_cPARP', 'Obs_MLKL'])
    td_PARP = ct.calculate_time_delay(PARP_MLKL_signals[:,0], sims.tspan)
    td_MLKL = ct.calculate_time_delay(PARP_MLKL_signals[:,1], sims.tspan)

    #-----------Time Delay vs. procaspase 8-----------------
    if (td_PARP is not None) & (td_MLKL is not None):
        delta_td.append(td_MLKL[0] - td_PARP[0])
        apopt_td.append(td_PARP[0])
        necro_td.append(td_MLKL[0])

#------------Plot Results--------------------------------
pl.ion()
pl.figure('Cell Death Signals')
pl.plot(range_XIAP, apopt_td)
pl.plot(range_XIAP, necro_td)
pl.xlabel('XIAP [molecules/cell]')
pl.ylabel('Time delay [s]')
pl.title('Time delay of apoptotic and necrotic signals vs. XIAP content')

pl.figure('Cell Death Signals 2')
pl.ion()
pl.plot(range_XIAP, delta_td)
pl.xlabel('XIAP  [molecules/cell]')
pl.ylabel('Difference in time delay between apoptotic and necrotic signals [s]')
pl.title('Time delay of apoptotic and necrotic signals vs. XIAPcontent')
pl.show()