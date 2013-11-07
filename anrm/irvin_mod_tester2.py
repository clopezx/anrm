"""
    Overview
    ========
    
    PySB implementations of the apoptosis-necrosis reaction model version 1.0
    (ANRM 1.0) originally published in [Irvin,NECRO2013]_.
    
    This model provides information about the dynamic biomolecular events that
    commit a cell to apoptosis or necrosis in response to TNFa or Fas signalling.
    PARP1 () serves, in this case, as a marker for apoptosis or necrosis. In
    apoptosis, PARP1 is cleaved whereas in necroptosis it is activated [].
    
    This file contains functions that implement the extrinsic apoptosis pathway
    and activation of a necroptosis reporter protein, PARP. 
    The modules are organized into three overall sections:
    - Receptor signalling and Bid activation.
    - Pore formation
    - PARP activitation/cleavage.

    These sections are further devided into...
    For receptor signalling and Bid activation:
    -- CD95 Ligation to formation of secondary complex
    -- TNFR1 ligation to formation of complex II
    -- Secondary complexes to formation of Riptosomes and Necrosomes and Bid
    activation.
    -- Bid translocation of the mitochondria and inhibition of anti-apoptotics
    For pore formation
    -- Lopez Pore formation (Adapted from Lopez Module)
    For PARP activation/cleavage
    -- Pore to PARP
    -- RIP1 to PARP
    """

import numpy
from pysb import *
from pysb.util import alias_model_components
from pysb.macros import *

#from shared_anrm import *
#from earm.shared import *

Model()

Parameter('KF', 1e-6) # Generic association rate constant
Parameter('KR', 1e-3) # Generic dessociation rate constant
Parameter('KC', 1)    # Generic catalytic rate constant
Parameter('KF2', 0)


def TNFR1_to_SecondaryComplex_monomers():
    """ Declares TNFa, TNFR1, TRADD, CompI, RIP1, A20, CYLD, NEMO and NFkB.
    Upon activation, TNFR1 gets endocytosed and post translationally modified After Complex I 
    has released TRADD and RIP1 it possibly gets recycled. This is represented by giving TNFR1 
    two states: norm and spent. TRADD has two states. CompI has two states. RIP1 has three states:
    Ub, PO4 and inactive. RIP1 binds FADD, Complex I, Bid-P and RIP3 (see SecondaryComplex_Bid).
    A20, CYLD and NEMO catalyze transformations of CompI and RIP1, maybe theycan be represented in the rates. 
    """
    
    Monomer('CompI', ['bDD', 'state'], {'state':['unmod', 'mod']}) #Neglecting RIP1:proC8 binding.. for simplicity.
    Monomer('RIP1', ['bDD', 'bRHIM', 'state'], {'state':['unmod', 'ub', 'po4', 'trunc', 'N']})

def TNFR1_to_SecondaryComplex():
    """Defines the interactoins from TNFR1 ligation to generation of secondary
    complexes as per ANRM 1.0.
        
    Uses TNFa, TNFR1, TRADD, CompI, RIP1 and NFkB. C8 active dimers and
    their associated parameters to generate rules that describe the ligand/receptor
    binding, FADD recruitment, proC8 and c-flip recruitment, activation of caspase
    and release of the secondary complex. This model assumes that one copy of proC8
    binds FADD before c-flip.
    
    This model converts proC8:proC8 to C8 (active caspase 8 dimer)
    This model also produces Secondary complex, FADD:proC8:c-Flip.

    RIP1-Ub recruitment to the CompI is not explicitly state in literature. But, this
    would be required to maintain RIP1 dependent TNFa signalling, if we allowed RIP1
    ubiquitination (inhibiting Apoptosis/Necrosis) to occur independently of TNFR1. 
    """
    
    Parameter('CompI_0' ,  1000) # complexes per cell
    Parameter('RIP1_0'  , 20000) # molecules per cell 20000

    Initial(CompI(bDD=None, state='mod'), CompI_0)              # Complex I
    Initial(RIP1(bDD=None, bRHIM = None, state = 'ub'), RIP1_0)   # RIP1

    # =========================================
    # TNFR1 ligation, formation of Complex I and release of RIP1 and TRADD rules
    # -----------------------------------------
    #   TNFa+ TNFR1 <-> TNFa:TNFR1
    #   TNFa:TNFR1 + TRADD <-> TNFa:TNFR1:TRADD >> CompI
    #   CompI + RIP1 <-> CompI:RIP1 >> [active]CompI:RIP1-Ub
    #   CompI + RIP1-Ub <-> CompI:RIP1-Ub 
    
    #   [active]CompI:RIP1-Ub >> NFkB # This reaction will consume the receptor.
    #   [active]CompI:RIP1-Ub >> [active]CompI:RIP1
    #   [active]CompI:RIP1 >> [active]CompI # A20 mediated degradation of RIP1
    #   [active]CompI:RIP1 >> [active]CompI + RIP1
    #   [active]CompI >> [active]TRADD + [norm]TNFR1 #receptor recycle typically distroys the ligand.

    #   RIP1 >> RIP1-Ub #These reactions were added because Doug Green reported that FADD and RIP1 bind
    #   RIP1-Ub >> RIP1 #independently of receptor, and FADD:RIP1 formation leads to Caspase 8 activation.
                        #These reaction decrease the amount of RIP1 by spontaneously ubiquitinating it. 

    # ------------------------------------------
    
    # -------------Complex I assembly----------------

    bind(CompI(bDD=None, state = 'mod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='ub'), 'bDD',[KF, KR])
    Rule('RIP1_Complex', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub') >> CompI(bDD=None, state='mod')%RIP1(bDD=ANY, bRHIM = None, state = 'unmod'),KR)
    Rule('RIP1_deg', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='unmod') >> CompI(bDD=None, state='mod'),KF)
    # Rule('RIP1_rel', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='unmod') >> CompI(bDD=None, state='mod') + RIP1(bDD=None, bRHIM = None,  state = 'unmod'), KF)
    bind(CompI(bDD=None, state='mod'), 'bDD', RIP1(bDD=None, bRHIM = None,  state = 'unmod'), 'bDD', [KF2, KR])

        
TNFR1_to_SecondaryComplex_monomers()
TNFR1_to_SecondaryComplex()


Observable('ComplexI', CompI(bDD=None))
Observable('ComplexI_RIP1mod', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub'))
Observable('ComplexI_RIP1unmod', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='unmod'))
Observable('Obs_RIP1', RIP1(bDD=None, state = 'unmod'))
Observable('Obs_RIP1_Ub', RIP1(bDD=None, state = 'ub'))
Observable('ComplexI_spent', CompI(bDD=None, state = 'unmod'))

