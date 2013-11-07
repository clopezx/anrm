"""
Apoptosis-Necrosis Reaction Network Model: 

"""

import numpy as np
import pylab as p

from pysb.integrate  import odesolve
#from anrm import irvin_AN_crosstalk
from anrm import irvin_anrm_experiment_8

m = irvin_anrm_experiment_8.model
#m = irvin_AN_crosstalk.model
t = np.linspace(0,10000,100)
yout = odesolve(m, t)

p.ion()
p.plot(t, yout['Obs_cPARP'], label = 'Cleaved Parp')
p.plot(t, yout['Obs_MLKL'], label = 'MLKL')

p.xlabel('time [sec]')
p.ylabel('PARP concentration [molecules per cell]')
p.legend()