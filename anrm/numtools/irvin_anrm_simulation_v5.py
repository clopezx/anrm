"""
Apoptosis-Necrosis Reaction Network Model: 

"""
import numpy as np
import pylab as p
import calibratortools as ct
import simulator_1_0 as sim
import pickle

from anrm.irvin_mod_v5_tester import model
from pysb.integrate  import odesolve

#-----------Calibrated Parameters-----------------------
#position = pickle.load(open('CompII_Hypthesis_123_addeddata_4run_v41_Position.pkl'))

#-----------Simulator Settings--------------------------
sims = sim.Settings()
sims.model = model
sims.tspan = np.linspace(0,20000,1000)
sims.estimate_params = model.parameters_rules()
sims.rtol = 1e-3
sims.atol = 1e-6

solve = sim.Solver(sims)
solve.run()

#-----------Initial Conditions--------------------------
ic_params  = model.parameters_initial_conditions()
conditions = ct.initial_conditions(['Bak_0', 'Bax_0', 'Bid_0', 'zVad_0'], [0.2e5, 40165, 0, 0], ic_params)
#20uM zVad == 9.6e6 zVad per cell for a cell volume of 8e-13m3
ysim = solve.simulate(position=None, observables=True, initial_conc = conditions)

yout = ct.extract_records(ysim, ['Obs_cPARP', 'Obs_MLKL','Obs_TNFa','Obs_NFkB', 'ComplexI','ComplexI_ub', 'ComplexI_TRAF', 'TRADD_RIP1', 'TRADD_RIP1_2','Obs_FADD_Sole', 'ComplexII','Bid_Trunc', 'Bid_PO4','Obs_RIP1', 'RIP1_Trunc', 'RIP3_Trunc', 'Necrosome','Obs_proC8', 'Obs_C8', 'Obs_C3ub', 'Obs_C3', 'Obs_pC3', 'RIP1_FADD','Obs_cPARP', 'Obs_PARP', 'Obs_MLKL','Obs_CytoC'])

p.ion()
p.plot(sims.tspan, yout[:,0], label = 'Cleaved Parp')
p.plot(sims.tspan, yout[:,1], label = 'MLKL')

p.xlabel('time [sec]')
p.ylabel('PARP and MLKL concentration [molecules per cell]')
p.legend()