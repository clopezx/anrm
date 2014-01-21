"""
Apoptosis-Necrosis Reaction Network Model: 

"""

import numpy as np
import pylab as p

import irvin_mod_v4_tester

from pysb.integrate  import odesolve

m = irvin_mod_v4_tester.model
t = np.linspace(0,10000,100)
yout = odesolve(m, t)

