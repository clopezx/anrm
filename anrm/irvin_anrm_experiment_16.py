from pysb import *
from pysb.util import alias_model_components
from earm import shared

from anrm import irvin_modules_v2 as irvin
from earm import lopez_modules as lopez 
from earm import albeck_modules as albeck
from anrm import merge

# -----Monomers-----
def compile_monomers():
    Model('m')
    # From irvin_modules
    irvin.Momomers_FasL_to_DISC()
    irvin.NFkB_cFlip_interaction_monomers()
    irvin.TNFR1_to_ComplexII_monomers()
    irvin.SecondaryComplex_to_Bid_monomers()
    irvin.rip1_to_MLKL_monmers()
    irvin.Momomers_zVad_to_C8()
    
    # From lopez_modules
    lopez.momp_monomers()

    # From albeck_modules
    albeck.apaf1_to_parp_monomers()
    return m.monomers

revised_monomers = {'Bid' :(['bf', 'state'], {'state':['U', 'T', 'M', 'po4']}),
    'C6'  :(['bf1','bf2', 'state'], {'state':['pro', 'A']}),
    'PARP':(['bf', 'state'], {'state':['U', 'C', 'A']}),
    'C3': (['bf', 'state'], {'state':['pro', 'A', 'ub', 'I']})}

monomer_edits = merge.Edit_Monomers(compile_monomers(), revised_monomers)
merged_monomers = monomer_edits.merged_monomers

Model('model')
model.monomers = merged_monomers
"""
    Apoptosis-Necrosis Reaction Network Model:
    
    Experimental conditions:
    1) TNFa sensitive cells.
    2) cFlip competes with procaspase for binding to FADD
    3) NFkB mediated cFlip expression
    4) cFlip-L can cleave Bid, but with lower efficiency
    5) RIP1 Hypothesis 2: RIP1-unmod =! RIP1-deub
        5a) RIP1 de-ubiquitination in the 2o complex
        5a) TRADD can bind RIP1-unmod as well as RIP1-deub
    6) Bid-po4 sequesters RIP1
    7) Z-VAD inhibits active C8 and C3 by irreversibly binding to their catalytic sites.
    
    """
# From irvin_modules
irvin.Initials_Fas_to_DISC()
irvin.TNFR1_to_ComplexII_Initials()
irvin.SecondaryComplex_to_Bid_initials()
irvin.rip1_to_MLKL_initials()
irvin.Initials_zVad_to_C8() #7

irvin.TNFR1_to_ComplexI() #1
irvin.FADD_to_C8()
irvin.cFlip_competitive_inhibition() #2
irvin.NFkB_cFlip_interaction() #3
irvin.cFLIP_L_Bid_interaction() #4
irvin.RIP1_Hypothesis_2()#5
irvin.ComplexI_to_NFkB()
irvin.ComplexII_Hypothesis_2() #5
irvin.TRADD_RIP1unmod_Hypothesis_2()#5a
irvin.RIP1_deubiqutination_Hypothesis_1()#5a
irvin.RIP1_to_SecondaryComplex()
irvin.RIP1_truncation()
irvin.Bid_Hypothesis() #6
irvin.C8_catalyzed_truncations()
irvin.rip1_to_MLKL()
irvin.zVad_to_C8() #7

# From lopez_modules
lopez.declare_initial_conditions()
lopez.translocate_tBid_Bax_BclxL()
lopez.tBid_activates_Bax_and_Bak()
lopez.tBid_binds_all_anti_apoptotics()
lopez.sensitizers_bind_anti_apoptotics()
lopez.effectors_bind_anti_apoptotics()
lopez.lopez_pore_formation(do_pore_transport=True)

# From albeck_modules
albeck.pore_to_parp()

# Observables
irvin.observables()


