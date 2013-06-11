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




# SECTION ONE: Receptor signalling and Bid Activation
# ===================================================
# This contains CD95_to_SecondaryComplex,
# TNFR1_to_SecondaryComplex,
# SecondaryComplex_to_Bid
# And the Shared Functions that describe Bid, Bad, Bax
# and Noxa binding to Anti-Apoptotic molecules.

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
    alias_model_components()

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

    Parameter('Fas_0'   ,   3000) # 3000 corresponds to 50ng/ml Fas
    Parameter('CD95_0'  ,    200) # 200 receptors per cell
    Parameter('FADD_0'  ,  1.0e3) # molecules per cell (arbitrarily assigned)1000
    Parameter('flip_L_0',  1.0e4) # molecules per cell
    Parameter('flip_S_0',  1.0e4) # molecules per cell
    Parameter('proC8_0' ,  2.0e4) # procaspase 8 molecules per cell 20000
    Parameter('C8_0'    ,      0) # active caspase 8 dimers per cell.
    Parameter('Bar_0'   ,  1.0e3) # Bar molecules per cell.
    alias_model_components()
    
    Initial(Fas(blig=None), Fas_0)       #Fas Ligand
    Initial(CD95(blig=None, bDD=None), CD95_0)     #Fas Receptor (CD95)
    Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0) #FADD
    Initial(flip_L(bDED=None), flip_L_0)   #c-Flip[L]
    Initial(flip_S(bDED=None), flip_S_0)   #c-Flip[S]
    Initial(proC8(bDED=None), proC8_0)    #procaspase 8
    Initial(C8(bC8=None), C8_0)       #caspase 8
    Initial(Bar(bC8=None), Bar_0)     #bifunctional apoptosis regulator
    
    # =========================================
    # CD95 ligation and formation of Secondary Complex rules
    # -----------------------------------------
    #   Fas + CD95 <-> Fas:CD95
    #   Fas:CD95 + FADD <-> Fas:CD95:FADD
    #   Fas:CD95:FADD + proC8 <-> Fas:CD95:FADD:proC8
    
    #   Fas:CD95:FADD:proC8 + proC8 <-> Fas:CD95:FADD:proC8:proC8 -> Fas:CD95:FADD + C8
    #   Fas:CD95:FADD:proC8 + flip_L <-> Fas:CD95:FADD:proC8:flip_L
    #   Fas:CD95:FADD:proC8 + flip_S <-> Fas:CD95:FADD:proC8:flip_S
    
    #   Fas:CD95:FADD:proC8:flip_L <-> Fas:CD95 + FADD:proC8:flip_L
    #   Fas:CD95:FADD:proC8:flip_S <-> Fas:CD95 + FADD:proC8:flip_S
    # ------------------------------------------
    
    # ==========================================
    # Regulation (downstream of CD95) of Caspase 8 activity
    # ------------------------------------------
    # C8 + BAR <-> C8:BAR
    # C8 + TRAF2 >> TRAF2
    
    # -------------DISC assembly----------------
    Parameter('kf1', 1e-6) # Generic association rate constant
    Parameter('kr1', 1e-3) # Generic dessociation rate constant
    Parameter('kf2', 1e-6) # Generic association rate constant
    Parameter('kr2', 1e-3) # Generic dessociation rate constant
    Parameter('kf3', 1e-6) # Generic association rate constant
    Parameter('kr3', 1e-3) # Generic dessociation rate constant
    Parameter('kf4', 1e-6) # Generic association rate constant
    Parameter('kr4', 1e-3) # Generic dessociation rate constant
    Parameter('kf5', 1e-6) # Generic association rate constant
    Parameter('kr5', 1e-3) # Generic dessociation rate constant
    Parameter('kf6', 1e-6) # Generic association rate constant
    Parameter('kr6', 1e-3) # Generic dessociation rate constant
    Parameter('kc6', 1e-1) # Generic catalytic rate constant
    Parameter('kf7', 1e-6) # Generic association rate constant
    Parameter('kr7', 1e-3) # Generic dessociation rate constant
    Parameter('kf8', 1e-10)# Low, non-zero affinity for forward reaction because this is documented as a release reaction.
    Parameter('kr8', 1e-3) # Generic dessociation rate constant
    
    alias_model_components()
    bind(Fas(blig=None), 'blig',  CD95(blig = None, bDD = None), 'blig', [kf1, kr1])
    bind(CD95(blig = ANY, bDD = None), 'bDD', FADD(bDD = None, bDED1 =None, bDED2 = None), 'bDD', [kf2, kr2])
    bind(FADD(bDD = ANY, bDED1 = None, bDED2 = None),'bDED1', proC8(bDED = None), 'bDED', [kf3, kr3])
    
    # For simplicity allow proC8 to bind FADD before any c-Flip do.
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', flip_L(bDED = None), 'bDED', [kf4, kr4])
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', flip_S(bDED = None), 'bDED', [kf5, kr5])
    
    
    # procaspase 8 dimerization and activation
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', proC8(bDED = None), 'bDED', [kf6, kr6])
    DISC_proC8 = CD95(blig=ANY, bDD=ANY) % Fas(blig=ANY) % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % proC8(bDED=ANY)%proC8(bDED=ANY)
    DISC = CD95(blig=ANY, bDD=ANY) % Fas(blig=ANY) % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)
    Rule('C8_activation', FADD(bDED1 = ANY, bDED2 = ANY)%proC8(bDED=ANY)%proC8(bDED = ANY) >> FADD(bDED1 = None, bDED2 = None) + C8(bC8 = None), kc6)

    # caspase 8 inhibition by BAR
    bind(Bar(bC8 = None), 'bC8', C8(bC8 = None), 'bC8', [kf7, kr7])
    
    # release of secondary complex from the DISC

    bind(FADD(bDD = None, bDED2 = ANY, bDED1 = ANY), 'bDD', CD95(blig = ANY, bDD=None), 'bDD', [kf8,kr8])

def TNFR1_to_SecondaryComplex_monomers():
    """ Declares TNFa, TNFR1, TRADD, CompI, RIP1, A20, CYLD, NEMO and NFkB.
    Upon activation, TNFR1 gets endocytosed and post translationally modified After Complex I
    has released TRADD and RIP1 it possibly gets recycled. This is represented by giving TNFR1
    two states: norm and spent. TRADD has two states. CompI has two states. RIP1 has three states:
    Ub, PO4 and inactive. RIP1 binds FADD, Complex I, Bid-P and RIP3 (see SecondaryComplex_Bid).
    A20, CYLD and NEMO catalyze transformations of CompI and RIP1, maybe be represented in the rates.
    """
    
    Monomer('TNFa', ['blig'])
    Monomer('TNFR1', ['blig', 'bDD', 'state'], {'state':['norm','spent']})
    Monomer('TRADD', ['bDD1', 'bDD2', 'state'], {'state':['active', 'inactive']})
    Monomer('CompI', ['bDD', 'state'], {'state':['unmod', 'mod']}) #Neglecting RIP1:proC8 binding.. for simplicity.
    Monomer('RIP1', ['bDD', 'bRHIM', 'bPARP', 'state'], {'state':['unmod', 'ub', 'po4', 'trunc', 'N']})
    Monomer('NFkB', ['bf'])
    alias_model_components()

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
    
    Parameter('TNFa_0'  ,  3000) # 3000 corresponds to 50ng/ml TNFa
    Parameter('TNFR1_0' ,   200) # 200 receptors per cell
    Parameter('TRADD_0' ,  1000) # molecules per cell (arbitrarily assigned)1000
    Parameter('RIP1_0'  , 20000) # molecules per cell 20000
    Parameter('CompI_0' ,     0) # complexes per cell
    Parameter('NFkB_0'  ,     0) # molecules per cell
    alias_model_components()
    
    Initial(TNFa(blig=None), TNFa_0)                                 # TNFa Ligand
    Initial(TNFR1(blig=None, bDD=None, state='norm'), TNFR1_0)       # TNFR1
    Initial(TRADD(bDD1=None, bDD2=None, state='inactive'), TRADD_0)  # TRADD
    Initial(CompI(bDD=None, state='unmod'), CompI_0)                 # Complex I
    Initial(RIP1(bDD=None, bRHIM = None, bPARP = None, state = 'ub'), RIP1_0)   # RIP1
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
    Parameter('kf9', 1e-6) # Generic association rate constant
    Parameter('kr9', 1e-3) # Generic dessociation rate constant
    Parameter('kf10', 1e-6) # Generic association rate constant
    Parameter('kr10', 1e-3) # Generic dessociation rate constant
    Parameter('kc10', 1)    # Generic catalytic rate constant
    Parameter('kf11', 1e-6) # Generic association rate constant
    Parameter('kr11', 1e-3) # Generic dessociation rate constant
    Parameter('kf12', 1e-6) # Generic association rate constant
    Parameter('kr12', 1e-3) # Generic dessociation rate constant
    Parameter('kc13', 1)    # Generic catalytic rate constant
    Parameter('kc14', 1)    # Generic catalytic rate constant
    Parameter('kc15', 1)    # Generic catalytic rate constant
    Parameter('kc16', 1)    # Generic catalytic rate constant
    Parameter('kr17', 1)    # Generic dissociation rate constant
    Parameter('kf17', 1e-10)# low affinity association rate constant because this reaction should net in release.
    Parameter('kc18', 1)    # Generic catalytic rate constant
    Parameter('kc19', 1e-5)    # Generic catalytic rate constant
    Parameter('kc20', 1e-7)    # Generic catalytic rate constant
    Parameter('kc21', 1e-10)    # Generic catalytic rate constant
    alias_model_components()

    bind(TNFa(blig=None), 'blig', TNFR1(blig=None, bDD=None, state='norm'), 'blig', [kf9, kr9])
    bind(TNFR1(blig = ANY, bDD = None, state =  'norm'), 'bDD', TRADD(bDD1=None, bDD2=None, state='inactive'), 'bDD1', [kf10, kr10])
    preCompI = TNFa(blig=ANY)%TNFR1(blig=ANY, bDD=ANY, state = 'norm')%TRADD(bDD1 = ANY, bDD2=None, state = 'inactive')
    Rule('CompI_formation', preCompI >> CompI(bDD=None, state = 'unmod'), kc10)
    
    # --------------Complex I - RIP1 Modification-----------------
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='unmod'), 'bDD',[kf11, kr11])
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='ub'), 'bDD',[kf12, kr12])
    
    
    Rule('CompI_Ub', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'unmod')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), kc13)
    
    Rule('CompI_Ub2', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), kc14)
    
    Rule('CompI_deUb', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub')>>CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,state='unmod'), kc15)

    Rule('RIP1_deg', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='unmod') >> CompI(bDD=None, state='mod'),kc16)
    
    bind(CompI(bDD=None, state='mod'), 'bDD', RIP1(bDD=None, bRHIM = None,  state = 'unmod'), 'bDD', [kf17, kr17])
    
    Rule('TNFR1_recycle', CompI(bDD=None, state='mod') >> TRADD(bDD1=None, bDD2 = None, state='active') + TNFR1(blig = None, bDD = None, state =  'norm'), kc18)
    
    Rule('NFkB_expression', CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub') + NFkB(bf=None), kc19)
    # --------------RIP1 Ubiquitination---------------------------
    Rule('RIP1_Ub', RIP1(bDD=None, bRHIM = None, state='unmod')>> RIP1(bDD=None, bRHIM = None, state='ub'), kc20)

    Rule('RIP1_deUb', RIP1(bDD=None, bRHIM = None, state='ub')>> RIP1(bDD=None, bRHIM = None, state='unmod'), kc21)

def SecondaryComplex_to_Bid_monomers():
    Monomer('Bid', ['bf', 'state'], {'state':['unmod', 'po4', 'trunc','M']})
    Monomer('BidK', ['bf']) #unknown Bid-kinase
    Monomer('RIP3', ['bRHIM', 'state'], {'state':['unmod', 'po4', 'trunc', 'N']})
    alias_model_components()

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
    Parameter('Bid_0'   , 2.0e4) # 2.0e4 molecules per cell
    Parameter('BidK_0'  , 5.0e3) # molecules per cell
    alias_model_components()
    
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
    Parameter('kf22', 1e-6) # Generic association rate constant
    Parameter('kr22', 1e-3) # Generic dessociation rate constant
    Parameter('Ka_RIP1_FADD',   1e-7) # Biochemica et Biophysica Acta 1834(2013) 292-300
    Parameter('Kd_RIP1_FADD',   1e-8) # Biochemica et Biophysica Acta 1834(2013) 292-300
    Parameter('kf24', 1e-6) # Generic association rate constant
    Parameter('kr24', 1e-3) # Generic dessociation rate constant
    Parameter('kc25', 1e-1)    # Generic catalytic rate constant
    Parameter('kf26', 1e-6) # Generic association rate constant
    Parameter('kr26', 1e-3) # Generic dessociation rate constant
    Parameter('kc26', 1e-1) # Generic catalytic rate constant
    Parameter('kc27', 1e-1) # Generic catalytic rate constant
    Parameter('kc28', 1e-1) # Generic catalytic rate constant
    Parameter('kf29', 1e-6) # Generic association rate constant
    Parameter('kr29', 1e-3) # Generic dessociation rate constant
    Parameter('kf30', 1e-6) # Generic association rate constant
    Parameter('kr30', 1e-3) # Generic dessociation rate constant
    Parameter('kf31', 1e-6) # Generic association rate constant
    Parameter('kr31', 1e-3) # Generic dessociation rate constant
    Parameter('kc31', 1e-1) # Generic catalytic rate constant
    Parameter('kf32', 1e-6) # Generic association rate constant
    Parameter('kr32', 1e-3) # Generic dessociation rate constant
    Parameter('kc32', 1e-1) # Generic catalytic rate constant
    Parameter('kf33', 1e-6) # Generic association rate constant
    Parameter('kr33', 1e-3) # Generic dessociation rate constant
    Parameter('kc33', 1e-1) # Generic catalytic rate constant
    Parameter('kf34', 1e-6) # Generic association rate constant
    Parameter('kr34', 1e-3) # Generic dessociation rate constant
    Parameter('kf35', 1e-6) # Generic association rate constant
    Parameter('kr35', 1e-3) # Generic dessociation rate constant
    Parameter('kc36', 1e-1) # Generic catalytic rate constant
    alias_model_components()
    
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', TRADD(bDD1=None, state = 'active'), 'bDD1', [kf22, kr22])
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'unmod'), 'bDD', [Ka_RIP1_FADD, Kd_RIP1_FADD])
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bDD', [kf24, kr24])
    # For simplicity, I am neglecting the binary intereaction that occurs between proC8 and RIP1.
    # Binding of proC8 and c-flip to FADD is accomplished in CD95_to_Secondary complex. 

    #--------------RIP1 Truncation reactions-------------
    #---Truncation by C8---------------------------------
    CIIA = TRADD(bDD2 = None, bDD1 = ANY, state = 'active') %  FADD(bDD=ANY, bDED1=None, bDED2=None)
    RIP_CIIA_proC8 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% TRADD(bDD2 = None, bDD1 = ANY, state = 'active') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    RIP_CIIB_proC8 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    Rule('RIP1_truncation_CIIA', RIP_CIIA_proC8 >> CIIA + C8(bC8=None) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), kc25)
    Rule('RIP1_truncation_CIIB', RIP_CIIB_proC8 >> FADD(bDD=None, bDED1=None, bDED2=None)+ C8(bC8=None) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), kc25)
    catalyze_state(C8(bC8=None), 'bC8', RIP1(bDD=None), 'bRHIM', 'state', 'unmod', 'trunc', [kf26, kr26, kc26])

    #---Truncation by proC8:cFlip_L---------------------
    Riptosome_FADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%FADD(bDD=1, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Rule('RIP1_truncation2', Riptosome_FADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), kc27)

    Riptosome_TRADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%TRADD(bDD1=ANY, bDD2=1)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Rule('RIP1_truncation1', Riptosome_TRADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), kc28)
    
    # -------------RIP3 Binding Interactions----------------
    Ripto1_Flip_S = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Necrosome1 = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=6, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 6, state = 'unmod')
    Rule('RIP3_binding1', Ripto1_Flip_S + RIP3(bRHIM= None, state = 'unmod') <> Necrosome1, kf29, kr29)
    
    Ripto2_Flip_S = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Necrosome2 = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=5, state='unmod') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 5, state = 'unmod')
    Rule('RIP3_binding2', Ripto2_Flip_S + RIP3(bRHIM= None, state = 'unmod') <> Necrosome2, kf30, kr30)
    
    #RIP3 Truncation
    catalyze_state(C8(bC8=None), 'bC8', RIP3(), 'bRHIM', 'state', 'unmod', 'trunc', [kf31, kr31, kc31])

    # Bid Phosphorylation and Truncation
    catalyze_state(BidK(), 'bf', Bid(), 'bf', 'state', 'unmod', 'po4', [kf32, kr32, kc32])
    catalyze_state(C8(bC8=None), 'bC8', Bid(), 'bf', 'state', 'unmod', 'trunc', [kf33, kr33, kc33])

    # Bid-PO4 competing with RIP1 for binding to Complex II
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', Bid(bf = None, state = 'po4'), 'bf', [kf34, kr34])
    # Bid-PO4 sequestering RIP1
    bind(RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bRHIM', Bid(bf = None, state = 'po4'), 'bf', [kf35, kr35])

    # RIP1 degradation
    degrade(RIP1(bDD=None, bRHIM = None, state = 'unmod'), kc36)

# Shared functions
# ================

# Monomer and initial condition declarations
# ------------------------------------------

def momp_monomers():
    alias_model_components()
    """Declare the monomers for the Bcl-2 family proteins, Cyto c, and Smac.

    'bf' is the site to be used for all binding reactions (with the
    exception of Bax and Bak, which have additional sites used for
    oligomerization).

    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization. Most Bcl-2 proteins have the potential for both cytoplasmic
    and mitochondrial localization, with the exceptions of Bak and Bcl-2,
    which are apparently constitutively mitochondrial.
    """

    # **Activators**.
    # Bid, states: Untruncated, Truncated, truncated and Mitochondrial
    # Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})
    # **Effectors**
    # Bax, states: Cytoplasmic, Mitochondrial, Active
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']})
    # Bak, states: inactive and Mitochondrial, Active (and mitochondrial)
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state':['M', 'A']})
    # **Anti-Apoptotics**
    Monomer('Bcl2', ['bf'])
    Monomer('BclxL', ['bf', 'state'], {'state':['C', 'M']})
    Monomer('Mcl1', ['bf', 'state'], {'state':['M', 'C']})
    # **Sensitizers**
    Monomer('Bad', ['bf', 'state'], {'state':['C', 'M']})
    Monomer('Noxa', ['bf', 'state'], {'state': ['C', 'M']})

    # **Cytochrome C and Smac**
    Monomer('CytoC', ['bf', 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', ['bf', 'state'], {'state':['M', 'C', 'A']})

def declare_initial_conditions():
    alias_model_components()
    """Declare initial conditions for Bcl-2 family proteins, Cyto c, and Smac.
    """
    #Parameter('Bid_0'   , 4.0e4) # Bid
    Parameter('BclxL_0' , 2.0e4) # cytosolic BclxL
    Parameter('Mcl1_0'  , 2.0e4) # Mitochondrial Mcl1
    Parameter('Bcl2_0'  , 2.0e4) # Mitochondrial Bcl2
    Parameter('Bad_0'   , 1.0e3) # Bad
    Parameter('Noxa_0'  , 1.0e3) # Noxa
    Parameter('CytoC_0' , 5.0e8) # 5.0e5 cytochrome c
    Parameter('Smac_0'  , 1.0e8) # 1.0e5 Smac
    Parameter('Bax_0'   , 0.8e7) # 0.8e5 Bax
    Parameter('Bak_0'   , 0.2e7) # 0.2e5 Bak

    alias_model_components()

    #Initial(Bid(bf=None, state='U'), Bid_0)
    Initial(Bad(bf=None, state='C'), Bad_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
    Initial(Bcl2(bf=None), Bcl2_0)
    Initial(BclxL (bf=None, state='C'), BclxL_0)
    Initial(Mcl1(bf=None, state='M'), Mcl1_0)
    Initial(Noxa(bf=None, state='C'), Noxa_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)

def translocate_tBid_Bax_BclxL():
    alias_model_components()
    Parameter('Kf_transloca',   1e-1) # from Lopez_modules...
    Parameter('Kr_transloca',   1e-3)
    """tBid, Bax and BclXL translocate to the mitochondrial membrane."""
    
    equilibrate(Bid(bf=None, state='trunc'), Bid(bf=None, state='M'), [Kf_transloca, Kr_transloca])

    free_Bax = Bax(bf=None, s1=None, s2=None) # Alias for readability
    equilibrate(free_Bax(state='C'), free_Bax(state='M'),
                [1e-1, 1e-3])

    equilibrate(BclxL(bf=None, state='C'), BclxL(bf=None, state='M'),
                [1e-1, 1e-3])

def tBid_activates_Bax_and_Bak():
    alias_model_components()
    """tBid activates Bax and Bak."""
    Parameter('kf37', 1e-6) # Generic association rate constant
    Parameter('kr37', 1e-3) # Generic dessociation rate constant
    Parameter('kc37', 1e-1) # Generic catalytic rate constant
    catalyze_state(Bid(state='M'), 'bf', Bax(state='M'), 'bf', 'state', 'M', 'A', [1e-6, 1e-3, 1e-1])

    Parameter('kf38', 1e-6) # Generic association rate constant
    Parameter('kr38', 1e-3) # Generic dessociation rate constant
    Parameter('kc38', 1e-1) # Generic catalytic rate constant
    catalyze_state(Bid(state='M'), 'bf', Bak(state='M'), 'bf', 'state', 'M', 'A', [1e-6, 1e-3, 1e-1])

N_A = 6.022e23
V   = 1.0e-12
def tBid_binds_all_anti_apoptotics():
    alias_model_components()
    """tBid binds and inhibits Bcl2, Mcl1, and Bcl-XL.

    The entries given in the `bind_table` are dissociation constants taken
    from Certo et al. (see ref). Dissociation constants in Certo et al.
    were published as nanomolar binding affinities; here they are converted
    into units of numbers of molecules by multiplying by `N_A` (Avogadro's
    number) and `V` (a default cell volume, specified in :doc:`shared`.

    The default forward rate represents diffusion limited association
    (1e6 Molar^-1 s^-1) and is converted into units of molec^-1 s^-1 by dividing
    by `N_A*V`.

    Certo, M., Del Gaizo Moore, V., Nishino, M., Wei, G., Korsmeyer, S.,
    Armstrong, S. A., & Letai, A. (2006). Mitochondria primed by death signals
    determine cellular addiction to antiapoptotic BCL-2 family members. Cancer
    Cell, 9(5), 351-365. `doi:10.1016/j.ccr.2006.03.027`
    """
    # Doug Green's "MODE 1" inhibition
    bind_table([[                        Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bid(state='M'),  66e-9*N_A*V,       12e-9*N_A*V,      10e-9*N_A*V]],
               'bf', 'bf', kf=1e6/(N_A*V))

def sensitizers_bind_anti_apoptotics():
    alias_model_components()
    """Binding of Bad and Noxa to Bcl2, Mcl1, and Bcl-XL.

    See comments on units for :py:func:`tBid_binds_all_anti_apoptotics`.
    """

    bind_table([[                        Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bad(state='M'),  11e-9*N_A*V,       10e-9*N_A*V,             None],
                [Noxa(state='M'),        None,              None,      19e-9*N_A*V]],
               'bf', 'bf', kf=1e-6)

def effectors_bind_anti_apoptotics():
    alias_model_components()
    """Binding of Bax and Bak to Bcl2, BclxL, and Mcl1.

    Affinities of Bak for Bcl-xL and Mcl-1 are taken from Willis et al.

    Preferential affinity of Bax for Bcl-2 and Bcl-xL were taken from Zhai et
    al.  Bax:Bcl2 and Bax:Bcl-xL affinities were given order of magnitude
    estimates of 10nM.

    See comments on units for :py:func:`tBid_binds_all_anti_apoptotics`.

    Willis, S. N., Chen, L., Dewson, G., Wei, A., Naik, E., Fletcher, J. I.,
    Adams, J. M., et al. (2005). Proapoptotic Bak is sequestered by Mcl-1 and
    Bcl-xL, but not Bcl-2, until displaced by BH3-only proteins. Genes &
    Development, 19(11), 1294-1305. `doi:10.1101/gad.1304105`

    Zhai, D., Jin, C., Huang, Z., Satterthwait, A. C., & Reed, J. C. (2008).
    Differential regulation of Bax and Bak by anti-apoptotic Bcl-2 family
    proteins Bcl-B and Mcl-1. The Journal of biological chemistry, 283(15),
    9580-9586.  `doi:10.1074/jbc.M708426200`
    """

    bind_table([[                            Bcl2,  BclxL(state='M'),         Mcl1],
                [Bax(state='A'), 10e-9*N_A*V,       10e-9*N_A*V,         None],
                [Bak(state='A'),        None,       50e-9*N_A*V,  10e-9*N_A*V]],
               'bf', 'bf', kf=1e6/(N_A*V))

def lopez_pore_formation(do_pore_transport=True):
    """ Pore formation and transport process used by all modules.
    """
    alias_model_components()

    # Rates
    pore_max_size = 4
    pore_rates = [[2.040816e-04,  # 1.0e-6/v**2
                   1e-3]] * (pore_max_size - 1)
    pore_transport_rates = [[2.857143e-5, 1e-3, 10]] # 2e-6 / v?

    # Pore formation by effectors
    assemble_pore_sequential(Bax(bf=None, state='A'), 's1', 's2', pore_max_size, pore_rates)
    assemble_pore_sequential(Bak(bf=None, state='A'), 's1', 's2', pore_max_size, pore_rates)

    # CytoC, Smac release
    if do_pore_transport:
        pore_transport(Bax(bf=None, state='A'),'s1','s2','bf', 4, 4, CytoC(state='M'),
                      'bf', CytoC(state='C'), pore_transport_rates)
        pore_transport(Bax(bf=None, state='A'), 's1','s2','bf', 4, 4, Smac(state='M'),
                      'bf', Smac(bf=None, state='C'), pore_transport_rates)
        pore_transport(Bak(bf=None, state='A'), 's1','s2','bf', 4, 4, CytoC(state='M'),
                      'bf', CytoC(bf=None, state='C'), pore_transport_rates)
        pore_transport(Bak(bf=None, state='A'), 's1','s2','bf', 4, 4, Smac(state='M'),
                      'bf', Smac(bf=None, state='C'), pore_transport_rates)

def apaf1_to_parp_monomers():
    alias_model_components()
    """ Declares CytochromeC, Smac, Apaf-1, the Apoptosome, Caspases 3, 6, 9,
    XIAP and PARP.

    The package variable 'bf' specifies the name of the site to be used
    for all binding reactions.

    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """

    # Cytochrome C
    Monomer('Apaf', ['bf', 'state'], {'state':['I', 'A']}) # Apaf-1
    Monomer('Apop', ['bf']) # Apoptosome (activated Apaf-1 + caspase 9)
    # Csp 3, states: pro, active, ubiquitinated
    Monomer('C3', ['bf', 'state'], {'state':['pro', 'A', 'ub']})
    # Caspase 6, states: pro-, Active
    Monomer('C6', ['bf1','bf2', 'state'], {'state':['pro', 'A']})
    Monomer('C9', ['bf']) # Caspase 9
    # PARP, states: Uncleaved, Cleaved
    Monomer('PARP', ['bf', 'state'], {'state':['U', 'C', 'A']})
    Monomer('XIAP', ['bf']) # X-linked Inhibitor of Apoptosis Protein
    
def pore_to_parp():
    """Defines what happens after the pore is activated and Cytochrome C and
    Smac are released.

    Uses CytoC, Smac, Apaf, Apop, C3, C6, C8, C9, PARP, XIAP monomers and their
    associated parameters to generate the rules that describe apoptosome
    formation, XIAP inhibition, activation of caspases (including
    caspase-6-mediated feedback), and cleavage of effector caspase substrates
    as specified in EARM 1.0.

    Declares initial conditions for CytoC, Smac, Apaf-1, Apoptosome, caspases
    3, 6, and 9, XIAP, and PARP.
    """
    Parameter('Kf_Apaf_acti',   5e-7) # from Albeck_modules.py
    Parameter('Kf_Apop_asse',   5e-8) # from Albeck_modules.py
    Parameter('Kf_C3_activa',   5e-9) # from Albeck_modules.py
    Parameter('Kf_Apop_inhi',   2e-6) # from Albeck_modules.py
    Parameter('Kf_Smac_inhi',   7e-6) # from Albeck_modules.py
    Parameter('Kf_C3_activ2',   1e-7) # from Albeck_modules.py
    Parameter('Kf_C3_ubiqui',   2e-8) # from Albeck_modules.py (Adjusted from 2e-6)
    Parameter('Kc_C3_ubiqui',   1e-1) # from Albeck_modules.py
    Parameter('Kr_PARP_clea',   1e-2) # from Albeck_modules.py
    Parameter('Kf_C8_activ2',   7e-6) # It is 3e-8 in the Albeck_modules.py which activate C8 without first dimerizing it.
    Parameter('Kr_C8_activ2',   1) # Since, I added a dimerization step I have to adjust this perameter as well.
    Parameter('Kc_C8_activ2',   1e-1)
    Parameter('Kf_Bax_activ',   1e-7) # from Albeck_modules.py

    # Declare initial conditions:

    Parameter('Apaf_0'  , 1.0e6) # 1.0e5 Apaf-1
    Parameter('C3_0'    , 1.0e4) # procaspase-3 (pro-C3)
    Parameter('C6_0'    , 1.0e4) # procaspase-6 (pro-C6)
    Parameter('C9_0'    , 1.0e6) # 1.0e5 procaspase-9 (pro-C9)
    Parameter('XIAP_0'  , 1.0e5) # X-linked inhibitor of apoptosis protein
    Parameter('PARP_0'  , 1.0e6) # C3* substrate

    alias_model_components()

    Initial(Apaf(bf=None, state='I'), Apaf_0)
    Initial(C3(bf=None, state='pro'), C3_0)
    Initial(C6(bf1=None, bf2 = None, state='pro'), C6_0)
    Initial(C9(bf=None), C9_0)
    Initial(PARP(bf=None, state='U'), PARP_0)
    Initial(XIAP(bf=None), XIAP_0)

    # CytoC and Smac activation after release
    # --------------------------------------
    Parameter('ktr39', 1e-6) # Generic association rate constant
    Parameter('ktr40', 1e-6) # Generic dessociation rate constant
    equilibrate(Smac(bf=None, state='C'), Smac(bf=None, state='A'), [1e-6, 1e-6])
    equilibrate(CytoC(bf=None, state='C'), CytoC(bf=None, state='A'),[1e-6, 1e-6])

    # Apoptosome formation
    # --------------------
    #   Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
    #   aApaf + pC9 <-->  Apop
    #   Apop + pC3 <-->  Apop:pC3 --> Apop + C3

    Parameter('kr41', 1e-3) # Generic dessociation rate constant
    Parameter('kc41', 1e-1) # Generic catalytic rate constant
    catalyze_state(CytoC(state = 'A'), 'bf', Apaf(), 'bf', 'state', 'I', 'A', [5e-7, 1e-3, 1e-1])
    
    Parameter('kr42', 1e-3) # Generic dessociation rate constant
    bind(Apaf(bf=None, state='A'),'bf', C9(bf=None), 'bf', [5e-8, kr42])
    
    Parameter('kc43', 1e-1) # Generic catalytic rate constant
    Rule('Apoptosome', Apaf(bf=ANY, state = 'A')%C9(bf=ANY)>>Apop(bf=None), kc43)
    #one_step_conv.. I could not find this macro in the tutorial. I think it's removed.
    #one_step_conv(Apaf(state='A'), C9(), Apop(bf=None), [Kf_Apop_asse, KR])
    Parameter('kr44', 1e-3) # Generic dessociation rate constant
    Parameter('kc44', 1e-1) # Generic catalytic rate constant
    catalyze_state(Apop(), 'bf', C3(), 'bf', 'state', 'pro', 'A', [Kf_C3_activa, kr44, kc44])

    # Apoptosome-related inhibitors
    # -----------------------------
    #   Apop + XIAP <-->  Apop:XIAP  
    #   cSmac + XIAP <-->  cSmac:XIAP  
    Parameter('kr45', 1e-3) # Generic dessociation rate constant
    bind(Apop(bf=None), 'bf', XIAP(bf=None), 'bf', [Kf_Apop_inhi, kr45])
    Parameter('kr46', 1e-3) # Generic dessociation rate constant
    bind(Smac(bf=None, state='A'), 'bf', XIAP(bf=None), 'bf', [Kf_Smac_inhi, kr46])

    # Caspase reactions
    # -----------------
    # Includes effectors, inhibitors, and feedback initiators:
    #
    #   pC3 + C8 <--> pC3:C8 --> C3 + C8 CSPS
    #   pC6 + C3 <--> pC6:C3 --> C6 + C3 CSPS
    #   XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U CSPS
    #   PARP + C3 <--> PARP:C3 --> CPARP + C3 CSPS
    #   pC8 + C6 <--> pC8:C6 --> C8 + C6 CSPS
    
    Parameter('kf47', 1e-7) # Generic association rate constant
    Parameter('kr47', 1e-3) # Generic dessociation rate constant
    Parameter('kc47', 1e-1) # Generic catalytic rate constant
    catalyze_state(C8(), 'bC8', C3(),'bf', 'state', 'pro','A', [kf47, kr47, kc47])
    
    Parameter('kf48', 1e-6) # Generic association rate constant
    Parameter('kr48', 1e-3) # Generic dessociation rate constant
    Parameter('kc48', 1e-1) # Generic catalytic rate constant
    catalyze_state(XIAP(), 'bf', C3(), 'bf', 'state', 'A', 'ub', [kf48, kr48, kc48])
    
    
    Parameter('kf49', 1e-6) # Generic association rate constant
    Parameter('kr49', 1e-3) # Generic dessociation rate constant
    Parameter('kc49', 1e-1) # Generic catalytic rate constant
    catalyze_state(C3(state='A'), 'bf', PARP(state='U'),'bf', 'state', 'U', 'C', [kf49, kr49, kc49])
    
    Parameter('kf50', 1e-6) # Generic association rate constant
    Parameter('kr50', 1e-3) # Generic dessociation rate constant
    Parameter('kc50', 1e-1) # Generic catalytic rate constant
    catalyze_state(C3(state='A'), 'bf', C6(bf2 = None), 'bf1', 'state', 'pro', 'A', [kf50, kr50, kc50])
    
    
    Parameter('kf51', 1e-6) # Generic association rate constant
    Parameter('kr51', 1e-3) # Generic dessociation rate constant
    bind(C6(bf1 = None, bf2 = None, state = 'A'), 'bf1', proC8(bDED = None), 'bDED', [kf51, kr51])
    
    Parameter('kf52', 1e-6) # Generic association rate constant
    Parameter('kr52', 1e-3) # Generic dessociation rate constant
    Parameter('kc52', 1e-1) # Generic catalytic rate constant
    bind(C6(bf1 = ANY, bf2 = None, state = 'A'), 'bf2', proC8(bDED = None), 'bDED', [kf52, kr52])
    Rule('C8_activation_byC6', C6(bf1 = ANY, bf2 = ANY, state = 'A')%proC8(bDED = ANY)%proC8(bDED = ANY) >> C8(bC8=None) + C6(bf1=None, bf2=None, state = 'A'), kc52)
         
    #catalyze(C6(state='A'), C8(state='pro'), C8(state='A'), [Kf_C8_activ2, KR, KC])

def rip1_to_parp():
    """Defines a connection between RIP1 and RIP3 and PARP activity observed in
    necroptosis. S. Joaun-Lanhouet et al. 2012 Obsered a requirement for RIP1
    RIP3 for PARP-1 activation. PARP-1 initiated TRAIL induced necroptosis in
    acidic extracellular conditions.

    Uses RIP1, RIP3 and PARP monomers and their associated parameters to generate
    the rules that describe RIP1, RIP3 phosphorylation, and PARP-1 activation.
    """

    Parameter('kc53', 1e-1) # Generic catalytic rate constant
    Parameter('kf54', 1e-6) # Generic association rate constant
    Parameter('kr54', 1e-3) # Generic dessociation rate constant
    Parameter('kc54', 1e-1) # Generic catalytic rate constant
    Parameter('kf55', 1e-6) # Generic association rate constant
    Parameter('kr55', 1e-3) # Generic dessociation rate constant
    Parameter('kc55', 1e-1) # Generic catalytic rate constant

    alias_model_components()
    Rule('Rip_PO4lation', RIP1(bRHIM=ANY, state = 'unmod')%RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'), kc53)
    catalyze_state(RIP1(state='po4'), 'bPARP', PARP(), 'bf', 'state', 'U', 'A', [kf54, kr54, kc54])
    catalyze_state(PARP(state='A'), 'bf', PARP(), 'bf', 'state', 'U', 'A', [kf55, kr55, kc55])
    #Parameter('Kc_PARPactiv',   1e-10) # This likely multistep process is modeled via a one-step
    # catalysis reaction with slow rate coefficient.
    #Parameter('Kc_PARPautoa',   4e-4)
    #Rule('PARP_activata', RIP1(state = 'po4') + PARP(bf=None, state='U') >> RIP1(state = 'po4') + PARP(bf=None, state='A'), Kc_PARPactiv)
    
    #Rule('PARP_autoacti', PARP(bf=None, state='A') + PARP(bf=None, state = 'U') >> PARP(bf=None, state='A') + PARP(bf=None, state='A'), Kc_PARPautoa)

def observables():
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
    Observable('Obs_TRADD', TRADD(state = 'inactive'))
    Observable('Obs_TRADDa', TRADD(bDD1 = None, state = 'active'))
    Observable('SecondaryComplex', FADD(bDD=None, bDED1 = ANY, bDED2 = ANY))
    Observable('Complex_IIA', TRADD(bDD1=ANY, bDD2=None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Riptosome1', RIP1(bDD = ANY, bRHIM = None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Riptosome2', RIP1(bDD = ANY, bRHIM = None)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Obs_RIP1_Ub', RIP1(state = 'ub'))
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
    Observable('Obs_C3ub', C3(state = 'ub'))
    Observable('Obs_C3', C3(state = 'A'))
    Observable('Obs_pC3', C3(state = 'pro'))
    Observable('Obs_Apaf', Apaf(state = 'A'))
    Observable('Obs_Apop', Apop())
    Observable('Obs_Cyc', CytoC(bf=None, state='C'))
    Observable('Obs_cSmac', Smac(state = 'C'))
    Observable('RIP1_nucl', RIP1(bDD = None, bRHIM=ANY, state = 'N')%RIP3(bRHIM=ANY, state = 'N'))
    Observable('Obs_cPARP', PARP(state='C'))
    Observable('Obs_aPARP', PARP(state='A'))
    Observable('Obs_PARP', PARP(state='U'))
    Observable('RIP1_TRADD', RIP1()%TRADD())