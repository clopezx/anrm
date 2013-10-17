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
Parameter('KC2', 10)  # Generic catalytic rate constant
Parameter('KC3', 1e-6)# Generic catalytic rate constant
Parameter('KE', 1e-4) # Generic gene expression rate constant

Parameter('Ka_RIP1_FADD',   1e-7) # Biochemica et Biophysica Acta 1834(2013) 292-300
Parameter('Kd_RIP1_FADD',   1e-8) # Biochemica et Biophysica Acta 1834(2013) 292-300

Parameter('Kf_Apaf_acti',   5e-7) # from Albeck_modules.py
Parameter('Kf_Apop_asse',   5e-8) # from Albeck_modules.py
Parameter('Kf_C3_activa',   5e-9) # from Albeck_modules.py
Parameter('Kf_Apop_inhi',   2e-6) # from Albeck_modules.py
Parameter('Kf_Smac_inhi',   7e-6) # from Albeck_modules.py
Parameter('Kf_C3_activ2',   1e-7) # from Albeck_modules.py
Parameter('Kf_C3_ubiqui',   2e-6) # from Albeck_modules.py
Parameter('Kc_C3_ubiqui',   1e-1) # from Albeck_modules.py
Parameter('Kr_PARP_clea',   1e-2) # from Albeck_modules.py
Parameter('Kf_C8_activ2',   3e-8) # from Albeck_modules.py
Parameter('Kf_Bax_activ',   1e-7) # from Albeck_modules.py

Parameter('Kf_transloca',   1e-1) # from Lopez_modules...
Parameter('Kr_transloca',   1e-3)

Parameter('Kc_PARPactiv',   1e-10) # This likely multistep process is modeled via a one-step
                                   # catalysis reaction with slow rate coefficient.

Parameter('Kdeg', 1e-7)

# SECTION ONE: Receptor signalling and Bid Activation
# ===================================================
# This contains CD95_to_SecondaryComplex,
# TNFR1_to_SecondaryComplex,
# SecondaryComplex_to_Bid
# And the Shared Functions that describe Bid, Bad, Bax
# and Noxa binding to Anti-Apoptotic molecules.

def TNFR1_to_SecondaryComplex_monomers():
    """ Declares TNFa, TNFR1, TRADD, CompI, RIP1, A20, CYLD, NEMO and NFkB.
    Upon activation, TNFR1 gets endocytosed and post translationally modified After Complex I 
    has released TRADD and RIP1 it possibly gets recycled. This is represented by giving TNFR1 
    two states: norm and spent. TRADD has two states. CompI has two states. RIP1 has three states:
    Ub, PO4 and inactive. RIP1 binds FADD, Complex I, Bid-P and RIP3 (see SecondaryComplex_Bid).
    A20, CYLD and NEMO catalyze transformations of CompI and RIP1, maybe theycan be represented in the rates. 
    """
    Monomer('TNFa', ['blig'])
    Monomer('TNFR1', ['blig', 'bDD', 'state'], {'state':['norm','spent']})
    Monomer('TRADD', ['bDD1', 'bDD2', 'state'], {'state':['active', 'inactive']})
    Monomer('CompI', ['bDD', 'state'], {'state':['unmod', 'mod']}) #Neglecting RIP1:proC8 binding.. for simplicity.
    Monomer('RIP1', ['bDD', 'bRHIM', 'state'], {'state':['unmod', 'ub', 'po4', 'trunc', 'N']})
    Monomer('NFkB', ['bf'])

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
    
    Parameter('TNFa_0'  ,  5000) # 3000 corresponds to 50ng/ml TNFa
    Parameter('TNFR1_0' ,   200) # 200 receptors per cell
    Parameter('TRADD_0' ,  1000) # molecules per cell (arbitrarily assigned)1000
    Parameter('CompI_0' ,     0) # complexes per cell
    Parameter('RIP1_0'  , 20000) # molecules per cell 20000
    Parameter('NFkB_0'  ,     0) # molecules per cell
    
    Initial(TNFa(blig=None), TNFa_0)                                 # TNFa Ligand
    Initial(TNFR1(blig=None, bDD=None, state='norm'), TNFR1_0)       # TNFR1
    Initial(TRADD(bDD1=None, bDD2=None, state='inactive'), TRADD_0)  # TRADD
    Initial(CompI(bDD=None, state='unmod'), CompI_0)                 # Complex I
    Initial(RIP1(bDD=None, bRHIM = None, state = 'ub'), RIP1_0)   # RIP1
    Initial(NFkB(bf=None), NFkB_0)

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
    bind(TNFa(blig=None), 'blig', TNFR1(blig=None, bDD=None, state='norm'), 'blig', [KF, KR])
    bind(TNFR1(blig = ANY, bDD = None, state =  'norm'), 'bDD', TRADD(bDD1=None, bDD2=None, state='inactive'), 'bDD1', [KF, KR])
    preCompI = TNFa(blig=ANY)%TNFR1(blig=ANY, bDD=ANY, state = 'norm')%TRADD(bDD1 = ANY, bDD2=None, state = 'inactive')
    Rule('CompI_formation', preCompI >> CompI(bDD=None, state = 'unmod'), KC)
    
    # --------------Complex I - RIP1 Modification-----------------
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='unmod'), 'bDD',[KF, KR])
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='ub'), 'bDD',[KF, KR]) 
    
    Rule('CompI_Ub', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'unmod')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), KC)
    Rule('CompI_Ub2', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), KC)
    Rule('CompI_deUb', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub')>>CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,state='unmod'),KC)
    bind(CompI(bDD=None, state='mod'), 'bDD', RIP1(bDD=None, bRHIM = None,  state = 'unmod'), 'bDD', [Parameter('k2', 0), KR])
    #Rule('RIP1_rel', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='unmod') >> CompI(bDD=None, state='mod') + RIP1(bDD=None, bRHIM = None,  state = 'unmod'), KC3)
    #bind(CompI(bDD=None, state = 'mod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='ub'), 'bDD',[KF, KR])
    #Rule('RIP1_deg', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='unmod') >> CompI(bDD=None, state='mod'),KC3)
    Rule('TNFR1_recycle', CompI(bDD=None, state='mod') >> TRADD(bDD1=None, bDD2 = None, state='active') + TNFR1(blig = None, bDD = None, state =  'norm'), KC)
    Rule('NFkB_expression', CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub') + NFkB(bf=None), KE)
    # --------------RIP1 Ubiquitination---------------------------
    Rule('RIP1_Ub', RIP1(bDD=None, bRHIM = None, state='unmod')>> RIP1(bDD=None, bRHIM = None, state='ub'), KC2)
    # Rule('RIP1_deUb', RIP1(bDD=None, bRHIM = None, state='ub')>> RIP1(bDD=None, bRHIM = None, state='unmod'), KR)


def SecondaryComplex_to_Bid_monomers():
    Monomer('Bid', ['bf', 'state'], {'state':['unmod', 'po4', 'trunc','M']})
    Monomer('BidK', ['bf']) #unknown Bid-kinase
    Monomer('RIP3', ['bRHIM', 'state'], {'state':['unmod', 'po4', 'trunc', 'N']})


def SecondaryComplex_to_Bid():
    """Defines the interactoins from TRADD RIP1 and FADD to generation of secondary
        complexes and truncation of Bid as per ANRM 1.0.
        
        Uses FADD, proC8, C8, flip_L, flip_S, TRADD, RIP1, RIP3 and Bid, and their
        associated parameters to generate rules that describe the FADD recruitment,
        proC8 and c-flip recruitment to the secondary complex, activation of caspase
        truncation of RIP1 and RIP3 and phosphorylation of Bid. This model assumes
        that one copy of proC8 binds FADD before c-flip. Among other things...
        
        This model converts proC8:proC8 to C8 (active caspase 8 dimer)
        This model also produces Secondary complex, FADD:proC8:c-Flip.
        """
    Parameter('RIP3_0'  , 2.0e4) # molecules per cell
    Parameter('Bid_0'   , 2.0e4) # molecules per cell
    Parameter('BidK_0'  , 5.0e3) # molecules per cell
    
    Initial(RIP3(bRHIM = None, state = 'unmod'), RIP3_0)   # RIP3
    Initial(Bid(bf = None, state = 'unmod'), Bid_0)        # Bid
    Initial(BidK(bf = None), BidK_0)
    # ==============================================================
    # Assembly of Complex II, Riptosome and Necrosome
    # --------------------------------------------------------------
    #   FADD + TRADD[active] <-> FADD:TRADD[active]
    #   FADD + RIP1 <-> FADD:RIP1
    #   TRADD + RIP1 <-> TRADD:RIP1
    
    #   CD95_to_secondary complex contains the rules for recruitment of proC8 to FADD.
    #       (RIP1 or TRADD):FADD + proC8 <-> (RIP1 or TRADD):FADD:proC8
    #       (RIP1 or TRADD):FADD:proC8 + proC8 <-> (RIP1 or TRADD):FADD:proC8:proC8
    #       (RIP1 or TRADD):FADD:proC8 + flip_L <-> (RIP1 or TRADD):FADD:proC8:flip_L
    #       (RIP1 or TRADD):FADD:proC8 + flip_S <-> (RIP1 or TRADD):proC8:flip_S
    
    #   RIP1%ProC8%ProC8(in a complex) >> RIP1[trunc] + C8 + (remains of the complex)
    #   RIP1%ProC8%cFlip[L](in a complex) >> RIP1[trunc] + remains of the complex)
    #   RIP1%cFlip[S](in a complex) + RIP3 >> RIP1:RIP3(in a complex, i.e. necrosome)
    
    #   RIP1 + C8 <-> RIP1:C8 >> RIP1[trunc] + C8
    #   RIP3 + C8 <-> RIP3:C8 >> RIP3[trunc] + C8
    #   Bid + C8 <-> Bid:C8 >> Bid[trunc] + C8
    
    # -------------Assembling Complex II-----------------
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', TRADD(bDD1=None, state = 'active'), 'bDD1', [KF, KR])
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'unmod'), 'bDD', [Ka_RIP1_FADD, Kd_RIP1_FADD])
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bDD', [KF, KR])
    # For simplicity, I am neglecting the binary intereaction that occurs between proC8 and RIP1.
    # Binding of proC8 and c-flip to FADD is accomplished in CD95_to_Secondary complex.
    
def CD95_to_SecondaryComplex_monomers():
    """ Declares Fas ligand, CD95, FADD, Flip_L, Flip_S procaspase8 and Caspase 8.
                
    'bf' is the site to be used for all binding reactions.
    
    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """
    Monomer('Fas', ['blig'])    #Fas Ligand
    Monomer('CD95', ['blig', 'bDD'])           #Fas Receptor (CD95)
    Monomer('FADD', ['bDD', 'bDED1','bDED2'])    #FADD
    Monomer('flip_L', ['bDED'])   #c-Flip[L] binds FADD at bca1 or bca2
    Monomer('flip_S', ['bDED'])   #c-Flip[S] binds FADD at bca1 or bca2
    Monomer('proC8', ['bDED'])    #procaspase 8 binds FADD at bca1 or bca2
    Monomer('C8', ['bC8'])        #active caspase 8
    Monomer('Bar', ['bC8'])       #bifunctional apoptosis regulator

def CD95_to_SecondaryComplex():
    """Defines the interactoins from CD95 ligation to generation of secondary
    complexes as per ANRM 1.0.
    
    Uses Fas, CD95, FADD, flip_L, flip_S, proC8 monomers and C8 active dimers, and
    their assicated parameters to generate rules that describe the ligand/receptor
    binding, FADD recruitment, proC8 and c-flip recruitment, activation of caspase
    and release of the secondary complex. This model assumes that one copy of proC8
    binds FADD before c-flip.
    
    This model converts proC8:proC8 to C8 (active caspase 8 dimer)
    This model also produces Secondary complex, FADD:proC8:c-Flip.
    """
    
    Parameter('Fas_0'   ,    000) # 3000 corresponds to 50ng/ml Fas
    Parameter('CD95_0'  ,    200) # 200 receptors per cell
    Parameter('FADD_0'  ,  1.0e3) # molecules per cell (arbitrarily assigned)1000
    Parameter('flip_L_0',  1.0e4) # molecules per cell
    Parameter('flip_S_0',  1.0e4) # molecules per cell
    Parameter('proC8_0' ,  2.0e4) # procaspase 8 molecules per cell 20000
    Parameter('C8_0'    ,      0) # active caspase 8 dimers per cell.
    Parameter('Bar_0'   ,  1.0e3) # Bar molecules per cell.
    
    Initial(Fas(blig=None), Fas_0)       #Fas Ligand
    Initial(CD95(blig=None, bDD=None), CD95_0)     #Fas Receptor (CD95)
    Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0) #FADD
    Initial(flip_L(bDED=None), flip_L_0)   #c-Flip[L]
    Initial(flip_S(bDED=None), flip_S_0)   #c-Flip[S]
    Initial(proC8(bDED=None), proC8_0)    #procaspase 8
    Initial(C8(bC8=None), C8_0)       #caspase 8
    Initial(Bar(bC8=None), Bar_0)     #bifunctional apoptosis regulator

    
CD95_to_SecondaryComplex_monomers()
CD95_to_SecondaryComplex()
TNFR1_to_SecondaryComplex_monomers()
TNFR1_to_SecondaryComplex()

SecondaryComplex_to_Bid_monomers()
SecondaryComplex_to_Bid()

Observable('Obs_TNFa', TNFa(blig =  None))
Observable('Obs_Fas', Fas(blig = None))
Observable('Obs_TNFR1', TNFR1(blig = None))
Observable('Obs_CD95', CD95(blig = None))
Observable('CD95_Fas', CD95(blig = ANY))
Observable('DISC', CD95(blig = ANY, bDD=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2 = ANY))
Observable('FADD_proC8_proC8', FADD(bDED1 = ANY, bDED2 = ANY)%proC8(bDED=ANY)%proC8(bDED = ANY))
Observable('TNFR1_TNF', TNFR1(blig=ANY,bDD = ANY))
Observable('ComplexI', CompI())
Observable('Obs_RIP1', RIP1(state = 'unmod'))
Observable('Obs_RIP1_Ub', RIP1(state = 'ub'))
Observable('Obs_TRADD', TRADD(state = 'inactive'))
Observable('Obs_TRADDa', TRADD(bDD1=None, state = 'active'))
Observable('Obs_FADD', FADD(bDD = None, bDED1 = None, bDED2 = None))
Observable('FADD_TRADD', FADD()%TRADD())
Observable('RIP1_TRADD', RIP1()%TRADD())
Observable('Obs_RIP1so', RIP1(bDD=None, bRHIM = None,  state = 'unmod'))

Observable('SecondaryComplex', FADD(bDD=None, bDED1 = ANY, bDED2 = ANY))
Observable('Complex_IIA', TRADD(bDD1=ANY, bDD2=None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Riptosome1', RIP1(bDD = ANY, bRHIM = None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Riptosome2', RIP1(bDD = ANY, bRHIM = None)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Obs_RIP1_cFlip', FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY))
Observable('CompI_RIP1', CompI(bDD=ANY))
Observable('CompI_mod', CompI(state='mod'))
Observable('RIP1_Bid', RIP1()%Bid())
Observable('Bid_Riptosome1', Bid(bf= ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Bid_Riptosome2', Bid(bf= ANY)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
Observable('Bid_Trunc', Bid(state='trunc'))
Observable('Bid_PO4', Bid(state='po4'))
Observable('RIP1_Trunc', RIP1(state='trunc'))
Observable('RIP3_Trunc', RIP3(state='trunc'))
Observable('Necrosome', RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'))
Observable('Obs_proC8', proC8())
Observable('Obs_C8', C8())
#Observable('Obs_C3', C3(state = 'A'))
#Observable('Obs_Apaf', Apaf(state = 'A'))
#Observable('Obs_Apop', Apop())
#Observable('Obs_Cyc', CytoC(bf=None, state='C'))
#Observable('Obs_Smac', Smac(state = 'C'))
#Observable('RIP1_nucl', RIP1(bDD = None, bRHIM=ANY, state = 'N')%RIP3(bRHIM=ANY, state = 'N'))
#Observable('Obs_cPARP', PARP(state='C'))
#Observable('Obs_aPARP', PARP(state='A'))
#Observable('Obs_PARP', PARP(state='U'))"""

