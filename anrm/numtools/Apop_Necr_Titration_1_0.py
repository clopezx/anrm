# Runs ANRM 1.0 (Irvin et. al 2013) under a range of pro-apoptotic and pro-necrotic caspase 8
# concentrations.

import numpy as np
import pylab as pl
import pickle
import calibratortools as ct
import simulator_1_0 as sim

# ----------Model and Initial Conditions----------------
from anrm.irvin_mod_v5_tester  import model

range_proC8 = np.linspace(10,50000,101) #starting at zero causes NaNs when you normalize the data.
range_cFlip = np.linspace(0,200000,201)
range_Bid   = np.linspace(0,24000,101)
range_RIP1  = np.linspace(10,35000, 101)
range_BidK  = np.linspace(0, 10000, 51)

#-----------Calibrated Parameters-----------------------
position = pickle.load(open('CompII_Hypthesis_123_newtopology_2run_v40_Position.pkl'))

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,72000,4000) #24hrs converted to seconds (4000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-5
sims.atol = 1e-5

solve = sim.Solver(sims)
solve.run()

delta_td = []
apopt_td = []
necro_td = []

condition_variable = 'RIP1_0'
graph_name = 'RIP1'
rangecv = range_RIP1

for i in rangecv:
    #-----------Initial Conditions--------------------------
    ic_params  = model.parameters_initial_conditions()
    conditions = ct.initial_conditions([condition_variable], [i], ic_params)
    ysim = solve.simulate(position = position, observables=True, initial_conc = conditions)

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
pl.plot(rangecv, apopt_td)
pl.plot(rangecv, necro_td)
pl.xlabel('%s [molecules/cell]' % graph_name)
pl.ylabel('Time delay [s]')
pl.title('Time delay of apoptotic and necrotic signals vs. %s content' % graph_name)

pl.figure('Cell Death Signals 2')
pl.ion()
pl.plot(rangecv, delta_td)
pl.xlabel('%s  [molecules/cell]' % graph_name)
pl.ylabel('Difference in time delay between apoptotic and necrotic signals [s]')
pl.title('Time delay of apoptotic and necrotic signals vs. %s content' % graph_name)
pl.show()