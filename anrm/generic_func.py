"""
Apoptosis-Necrosis Reaction Network Model: 

"""
from pysb import *
from pysb.macros import *

Model()
Parameter('ka', 1e-6)
Parameter('kd', 1e-3)

Monomer('A', ['bf'])
Monomer('B', ['bf'])

Parameter('A_0', 100)
Parameter('B_0', 100)

Initial(A(bf=None), A_0)
Initial(B(bf=None), B_0)

bind(A(bf = None), 'bf', B(bf = None), 'bf', [ka, kd])

