# Runs ANRM 1.0 (Irvin et. al 2013) under a range of pro-apoptotic and pro-necrotic caspase 8
# concentrations.

import numpy as np
import pylab as p
import pickle
import calibratortools as ct
import simulator_1_0 as sim

# ----------Model and Initial Conditions----------------
from anrm.irvin_mod_v4_tester  import model

range_proC8 = [10, 4875, 9750, 19500, 39022, 78000] #starting at zero causes NaNs when you normalize the data.
range_cFlip = [0, 4875, 9750, 19500, 39022, 78000]
range_Bid   = [0, 4000, 8000, 12044, 16000, 20000]
range_RIP1  = [0, 4000, 8000, 12044, 16000, 20000]
range_BidK  = [0, 1000, 2500, 5000, 7500, 10000]


#-----------Calibrated Parameters-----------------------
position = pickle.load(open('CompII_Hypthesis_123_addeddata_4run_v41_Position.pkl'))

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,20000,2000) #10hrs converted to seconds (1000 timepoints)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-5
sims.atol = 1e-5

solve = sim.Solver(sims)
solve.run()

delta_td = []
apopt_td = []
necro_td = []

p.ion()
yout = []

condition_variable = 'Bid_0'
observable = 'Obs_cPARP'
graph_name = 'Bid'
rangecv = range_Bid

for i in rangecv:
    #-----------Initial Conditions--------------------------
    ic_params  = model.parameters_initial_conditions()
    conditions = ct.initial_conditions([condition_variable], [i], ic_params)
    ysim = solve.simulate(position = position, observables=True, initial_conc = conditions)
    
    #-----------Plot Parp and MLKL--------------------------
    yout.append(ct.extract_records(ysim, [observable]))

for j in range(len(rangecv)):
    p.plot(sims.tspan, yout[j], label = '%s %s per cell' % (rangecv[j], graph_name))

p.title('PARP cleavage in cells with varying initial %s concentrations' % graph_name)
p.xlabel('time [sec]')
p.ylabel('PARP concentration [molecules per cell]')
p.legend()