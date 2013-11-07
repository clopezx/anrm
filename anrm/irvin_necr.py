"""
Apoptosis-Necrosis Reaction Network Model: 

"""

from pysb import *
from earm import shared
from anrm import irvin_modules as irvin
from earm import lopez_modules as lopez 
from earm import albeck_modules as albeck


Model()

# -----Monomers-----
# From irvin_modules
irvin.CD95_to_SecondaryComplex_monomers()
irvin.TNFR1_to_SecondaryComplex_monomers()
irvin.SecondaryComplex_to_Bid_monomers()
"""alias_model_components() was placed in on the line following declaration of a parameter. This function exports the Parameter as a global variable. And it appears to be required before any rules can use the Parameters
    Line 69: Monomer('C8', ['bC8']) is changed to
    Monomer('C8', ['bf', 'state'], {'state':['A']}) 
    Line 104: Initial(C8(bC8=None), C8_0) is changed to
    Initial(C8(bf=None, state = 'A'), C8_0) 
    Line 161: Rule('C8_activation' ...+ C8(bC8 = None...) is changed to 
    Rule ('C8_activation' ...+C8(bf = None, state = 'A')
    Line 164: bind(Bar(bC8 = None), 'bC8', C8(bC8 = None), 'bC8', [kf7, kr7]) 
    is changed to bind(Bar(bC8 = None), 'bC8', C8(bf = None, state = 'A'), 'bf', [kf7, kr7])
    Line 293: Deleted. Monomer('Bid', ['bf', 'state'], {'state':['unmod', 'po4', 'trunc', 'M']})
    Line 310: Deleted. Parameter('Bid_0'   , 2.0e4)
    Line 315: Deleted. Initial(Bid(bf = None, state = 'unmod'), Bid_0) 
    Line 382: Rule('RIP1_truncation_CIIA', RIP_CIIA_proC8 >> CIIA + C8(bC8 = None) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), kc25) is changed to
        Rule('RIP1_truncation_CIIA', RIP_CIIA_proC8 >> CIIA + C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), kc25)
    Line 383: Rule('RIP1_truncation_CIIB', RIP_CIIB_proC8 >> FADD(bDD=None, bDED1=None, bDED2=None)+ C8(bC8=None) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), kc25) is changed to
    Line 384: catalyze_state(C8(bC8=None), 'bC8', RIP1(bDD=None), 'bRHIM', 'state', 'unmod', 'trunc', [kf26, kr26, kc26])
    is changed to catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP1(bDD=None), 'bRHIM', 'state', 'unmod', 'trunc', [kf26, kr26, kc26])
    Line 403:catalyze_state(C8(bC8=None), 'bC8', RIP3(), 'bRHIM', 'state', 'unmod', 'trunc', [kf31, kr31, kc31])
        is changed to catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP3(), 'bRHIM', 'state', 'unmod', 'trunc', [kf31, kr31, kc31])
    Line 406: catalyze_state(BidK(), 'bf', Bid(), 'bf', 'state', 'unmod', 'po4', [kf32, kr32, kc32]) is changed to
        catalyze_state(BidK(), 'bf', Bid(), 'bf', 'state', 'U', 'po4', [kf32, kr32, kc32])
    Line 407: catalyze_state(C8(bC8=None), 'bC8', Bid(), 'bf', 'state', 'unmod', 'trunc', [kf33, kr33, kc33])
    is changed to catalyze_state(C8(bf = None, state = 'A'), 'bf', Bid(), 'bf', 'state', 'unmod', 'trunc', [kf33, kr33, kc33])
    Line 407: catalyze_state(C8(bf = None, state = 'A'), 'bf', Bid(), 'bf', 'state', 'unmod', 'trunc', [kf33, kr33, kc33])
        is changed to
        catalyze_state(C8(bf = None, state = 'A'), 'bf', Bid(), 'bf', 'state', 'U', 'T', [kf33, kr33, kc33])
    """

# From lopez_modules
lopez.momp_monomers()
""" Line 99: Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']}) 
    is changed to Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M', 'po4']})
    """

# From albeck_modules
albeck.apaf1_to_parp_monomers()
""" Line 118: Monomer('C6', ['bf', 'state'], {'state':['pro', 'A']})
    is changed to Monomer('C6', ['bf1','bf2', 'state'], {'state':['pro', 'A']})
    Line 124: Monomer('PARP', ['bf', 'state'], {'state':['U', 'C']})
    is chanbed to Monomer('PARP', ['bf', 'state'], {'state':['U', 'C', 'A']})"""

# Rules
irvin.CD95_to_SecondaryComplex()
irvin.TNFR1_to_SecondaryComplex()
irvin.SecondaryComplex_to_Bid()

lopez.declare_initial_conditions()
lopez.translocate_tBid_Bax_BclxL()
lopez.tBid_activates_Bax_and_Bak()
lopez.tBid_binds_all_anti_apoptotics()
lopez.sensitizers_bind_anti_apoptotics()
lopez.effectors_bind_anti_apoptotics()
lopez.lopez_pore_formation(do_pore_transport=True)

albeck.pore_to_parp()
""" line 232: C6(bf=None...) changed to C6(bf1=None, bf2=None...)
    line 273: inserted. from pysb.macros import catalyze_state
    line 277: catalyze is changed to catalyze_state because C6 has 'bf1' instead of
    binding site 'bf'. 
    line 286: inserted. from pysb.macros import bind as bind2 earm.shared contains
    a bind function
    line 287and 288 replaced activation of C8 with: 
    bind2(C6(bf1 = None, bf2 = None, state = 'A'), 'bf1', proC8(bDED = None), 'bDED', [kf51, kr51])
    bind2(C6(bf1 = ANY, bf2 = None, state = 'A'), 'bf2', proC8(bDED = None), 'bDED', [kf52, kr52])
    Rule('C8_activation_byC6', C6(bf1 = ANY, bf2 = ANY, state = 'A')%proC8(bDED = ANY)%proC8(bDED = ANY) >> C8(bf = None, state = 'A') + C6(bf1=None, bf2=None, state = 'A'), kc52)"""

irvin.rip1_to_parp()

# Observables
irvin.observables()


    

