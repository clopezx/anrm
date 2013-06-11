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
"""alias_model_components() was placed in on the line following declaration of a
    parameter. This function exports the Parameter as a global variable. And it
    appears to be required before any rules can use the Parameters"""
irvin.CD95_to_SecondaryComplex_monomers()
irvin.TNFR1_to_SecondaryComplex_monomers()
irvin.SecondaryComplex_to_Bid_monomers()

# From lopez_modules
lopez.momp_monomers()
""" irvin.SecondaryComplex_to_Bid_monomers() has declares the monomer "Bid" as:
    Monomer('Bid', ['bf', 'state'], {'state':['unmod', 'po4', 'trunc','M']})

    lopez.momp_monomers() declares the monomer "Bid" as:
    Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})
    
    In accordance with Dr. Zinkel's hypothesis, Irvin has added a Bid phosphorylation 
    event. In order for use, lopez.momp_monomers was modified (i.e. deleted Bid Monomer line)
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
""" line 140: line deleted. Bid_0 was defined in irvin.SecondaryComplex_to_Bid()
    line 153: line deleted. Bid initial condition defined in irvin.Secondary... 
    line 169: Bid(... state= 'T') changed to Bid(... state = trunc)"""

albeck.pore_to_parp()
""" line 232: C6(bf=None...) changed to C6(bf1=None, bf2=None...)
    line 273: inserted. from pysb.macros import catalyze_state
    line 274: C8(state='A') changed to C8(). Since proC8 dimerizes and converts to
    caspase 8. The active state is implied.
    line 274: catalyze(... is changed to catalyze_state(...). The catalyze macro
    does not  allow you to specify a binding site. Instead, it assumes the binding
    site is 'bf'. The binding site for C8 is 'bC8'.
    line 277: catalyze is changed to catalyze_state because C6 has 'bf1' instead of
    binding site 'bf'. 
    line 286: inserted. from pysb.macros import bind as bind2 earm.shared contains
    a bind function"""

irvin.rip1_to_parp()

# Observables
irvin.observables()


    

