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
    
    These modules are to be run in cooperation with Lopez Momp_monomers and Albeck Apaf_to_Parp modules.
    
    Parameters and Initial Conditions were taken from Jurket Cell studies:
    Hua, F., M. Cornejo, M. Cardone, C. Stokes, D. Lauffenburger
    J Immunol 2005; 175:985-995
    --concentrations listed in this were converted to copies per cell using a cell
    volume of 8e-13L http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=104668&ver=2
    --
    """

import numpy
from pysb import *
from pysb.util import alias_model_components
from pysb.macros import *


# SECTION ONE: Receptor signalling and Bid Activation
# ===================================================
# This contains CD95_to_SecondaryComplex,
# TNFR1_to_SecondaryComplex,
# SecondaryComplex_to_Bid

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
    Monomer('flip_L', ['bDED'])   #c-Flip[L] binds FADD at DED2
    Monomer('flip_S', ['bDED'])   #c-Flip[S] binds FADD at DED2
    Monomer('proC8', ['bDED'])    #procaspase 8 binds FADD at DED1 or DED2
    Monomer('C8', ['bf', 'state'], {'state':['A']})        #active caspase 8
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

    Parameter('Fas_0'   ,   6000) # 6000 corresponds to 100ng/ml Fas (J Immunol 2005)
    Parameter('CD95_0'  ,   4817) # 4817receptors per cell (J Immunol 2005)
    Parameter('FADD_0'  ,   8030) # 8300 molecules per cell (J Immunol 2005)
    Parameter('flip_L_0',  390220) # 39022 molecules per cell (J Immunol 2005)
    Parameter('flip_S_0',  390220) # molecules per cell assummed both isoforms w/conc.
    Parameter('proC8_0' ,  16057) # procaspase 8 molecules per cell 16057 (J Immunol 2005)
    Parameter('C8_0'    ,      0) # active caspase 8 dimers per cell.
    Parameter('Bar_0'   ,  1.0e3) # Bar molecules per cell. (Hela Cells)
    alias_model_components()
    
    Initial(Fas(blig=None), Fas_0)       #Fas Ligand
    Initial(CD95(blig=None, bDD=None), CD95_0)     #Fas Receptor (CD95)
    Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0) #FADD
    Initial(flip_L(bDED=None), flip_L_0)   #c-Flip[L]
    Initial(flip_S(bDED=None), flip_S_0)   #c-Flip[S]
    Initial(proC8(bDED=None), proC8_0)    #procaspase 8
    Initial(C8(bf=None, state = 'A'), C8_0)       #caspase 8
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
    alias_model_components()
    bind(Fas(blig=None), 'blig',  CD95(blig = None, bDD = None), 'blig', [1.89e-7, 1e-4])
    bind(CD95(blig = ANY, bDD = None), 'bDD', FADD(bDD = None, bDED1 =None, bDED2 = None), 'bDD', [1.04e-6, 0.2]) # (J Immunol 2005)
    bind(FADD(bDD = ANY, bDED1 = None, bDED2 = None),'bDED1', proC8(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    
    # For simplicity allow proC8 to bind FADD before any c-Flip do.
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', flip_L(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', flip_S(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    
    
    # procaspase 8 dimerization and activation
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', proC8(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    
    Rule('C8_activation', FADD(bDED1 = ANY, bDED2 = ANY)%proC8(bDED=ANY)%proC8(bDED = ANY) >> FADD(bDED1 = None, bDED2 = None) + C8(bf = None, state = 'A'), Parameter('k1', 0.3)) # (J Immunol 2005)

    # caspase 8 inhibition by BAR
    bind(Bar(bC8 = None), 'bC8', C8(bf = None, state = 'A'), 'bf', [1e-6, 1e-3])
    
    # release of secondary complex from the DISC
    bind(FADD(bDD = None, bDED2 = ANY, bDED1 = ANY), 'bDD', CD95(blig = ANY, bDD=None), 'bDD', [1e-10, 1e-3])

def cFlip_competitive_inhibition():
    """Competitative inhibition of C8 activation. CD95_to_SecondaryComplex(): has inhibition by 1) decreasing ability to cleave Bid (i.e. cFlip_L) or 2) preventing activation (i.e. cFlip_S). This also added inhibition by preventing binding of proc8 and FADD. This is required to demonstrate increased inhibition of apoptosis at high concentrations."""
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = None), 'bDED1', flip_S(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = None), 'bDED1', flip_L(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)

def NFkB_cFlip_interaction_monomers():
    Monomer('Flip_degradase', ['bf'])
    alias_model_components()

def NFkB_cFlip_interaction():
    """Several papers mention NFkB mediated expression of c-Flip as an important mediator of apoptosis and necrosis."""
    Parameter('Flip_degradase_0', 0)
    alias_model_components()
    
    Initial(Flip_degradase(bf=None), Flip_degradase_0)
    
    Rule('NFkB_cFlipL', NFkB() >> NFkB() + flip_L(bDED=None), Parameter('NFkB_FlipL', 1e-2))
    Rule('NFkB_cFlipS', NFkB() >> NFkB() + flip_S(bDED=None), Parameter('NFkB_FlipS', 1e-2))
    
    Rule('NFkB_degradase', NFkB() >> NFkB() + Flip_degradase(bf=None), Parameter('Deg_flip', 1e-6))
    Rule('Deg_cFlipL', Flip_degradase(bf=None) + flip_L(bDED=None) >> Flip_degradase(bf=None), Parameter('deg_FlipL', 5e-6))
    Rule('Deg_cFlipS', Flip_degradase(bf=None) + flip_S(bDED=None) >> Flip_degradase(bf=None), Parameter('deg_FlipS', 5e-6))

def cFLIP_L_Bid_interaction():
    Rule('Bid_to_SecComp', FADD(bDD = None, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED = ANY)%proC8(bDED = ANY) + Bid(bf=None, state='U') <> FADD(bDD = 1, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED = ANY)%proC8(bDED = ANY)%Bid(bf=1, state='U'), Parameter('Bid_Flip_Bindf', 1e-6),Parameter('Bid_Flip_Bindr', 1e-1))
    
    Rule('tBid_SlowFormation', FADD(bDD = ANY, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED = ANY)%proC8(bDED = ANY)%Bid(bf=ANY, state='U') >> FADD(bDD = None, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED = ANY)%proC8(bDED = ANY) + Bid(bf=None, state='T'), Parameter('Bid_Flip_cat', 1e-1))

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
    Monomer('RIP1', ['bDD', 'bRHIM', 'bPARP', 'state'], {'state':['unmod', 'ub', 'po4', 'trunc', 'deub']})
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
    
    Parameter('TNFa_0'  ,  6000) # 6000 corresponds to 100ng/ml TNFa
    Parameter('TNFR1_0' ,  4800) # 4800 receptors per cell
    Parameter('TRADD_0' ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('RIP1_0'  , 12044) # molecules per cell (arbitrarily assigned) 12044
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
    bind(TNFa(blig=None), 'blig', TNFR1(blig=None, bDD=None, state='norm'), 'blig', [2e-7, 1e-4])
    bind(TNFR1(blig = ANY, bDD = None, state =  'norm'), 'bDD', TRADD(bDD1=None, bDD2=None, state='inactive'), 'bDD1', [1e-6, 0.2])
    preCompI = TNFa(blig=ANY)%TNFR1(blig=ANY, bDD=ANY, state = 'norm')%TRADD(bDD1 = ANY, bDD2=None, state = 'inactive')
    Rule('CompI_formation', preCompI >> CompI(bDD=None, state = 'unmod'), Parameter('k2', 1))
    
    # --------------Complex I - RIP1 Modification-----------------
    # Complex I recruits RIP1
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='unmod'), 'bDD',[1e-6, 1e-3])
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='ub'), 'bDD',[1e-6, 1e-3])
    
    # Complex I undergos modifications (e.g. Ubiquitination) and elicits NFkB signaling
    Rule('CompI_Ub', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'unmod')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), Parameter('k3', 1))
    Rule('CompI_Ub2', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), Parameter('k4', 1))
    Rule('NFkB_signaling', CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub') + NFkB(bf=None), Parameter('k5', 1e-3))

    # De-ubiquitination and degradation or release of RIP1 from complex I
    Rule('CompI_deUb', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub')>>CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,state='unmod'), Parameter('k6',1))
    Rule('RIP1_deg', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub') >> CompI(bDD=None, state='mod'), Parameter('k7',1))
    
    bind(CompI(bDD=None, state='mod'), 'bDD', RIP1(bDD=None, bRHIM = None,  state = 'unmod'), 'bDD', [1e-10,1])
        
    # Receptor recycle
    Rule('TNFR1_recycle', CompI(bDD=None, state='mod') >> TRADD(bDD1=None, bDD2 = None, state='active') + TNFR1(blig = None, bDD = None, state =  'norm'), Parameter('k8',1))
    
    # --------------RIP1 Ubiquitination---------------------------
    Rule('RIP1_Ub', RIP1(bDD=None, bRHIM = None, state='unmod')>> RIP1(bDD=None, bRHIM = None, state='ub'), Parameter('k9', 1e-7))

def TNFR1_to_SecondaryComplex_Alternate():
    """Defines the interactions from TNFR1 ligation to generation of secondary
        complexes as per ANRM 1.0.
        
        Uses TNFa, TNFR1, TRADD, CompI, RIP1 and NFkB. C8 active dimers and
        their associated parameters to generate rules that describe the ligand/receptor
        binding, FADD recruitment, proC8 and c-flip recruitment, activation of caspase
        and release of the secondary complex. This model assumes that one copy of proC8
        binds FADD before c-flip.
        
        This model converts proC8:proC8 to C8 (active caspase 8 dimer)
        This model also produces Secondary complex, FADD:proC8:c-Flip.
        
        This
        RIP1-Ub recruitment to the CompI is not explicitly state in literature. But, this
        would be required to maintain RIP1 dependent TNFa signalling, if we allowed RIP1
        ubiquitination (inhibiting Apoptosis/Necrosis) to occur independently of TNFR1.
        """
    
    Parameter('TNFa_0'  ,  000) # 6000 corresponds to 100ng/ml TNFa
    Parameter('TNFR1_0' ,  4800) # 4800 receptors per cell
    Parameter('TRADD_0' ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('RIP1_0'  , 12044) # molecules per cell (arbitrarily assigned) 12044
    Parameter('CompI_0' ,     0) # complexes per cell
    Parameter('NFkB_0'  ,     0) # molecules per cell
    alias_model_components()
    
    Initial(TNFa(blig=None), TNFa_0)                                 # TNFa Ligand
    Initial(TNFR1(blig=None, bDD=None, state='norm'), TNFR1_0)       # TNFR1
    Initial(TRADD(bDD1=None, bDD2=None, state='inactive'), TRADD_0)  # TRADD
    Initial(CompI(bDD=None, state='unmod'), CompI_0)                 # Complex I
    Initial(RIP1(bDD=None, bRHIM = None, bPARP = None, state = 'unmod'), RIP1_0)   # RIP1
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
    bind(TNFa(blig=None), 'blig', TNFR1(blig=None, bDD=None, state='norm'), 'blig', [2e-7, 1e-4])
    bind(TNFR1(blig = ANY, bDD = None, state =  'norm'), 'bDD', TRADD(bDD1=None, bDD2=None, state='inactive'), 'bDD1', [1e-6, 0.2])
    preCompI = TNFa(blig=ANY)%TNFR1(blig=ANY, bDD=ANY, state = 'norm')%TRADD(bDD1 = ANY, bDD2=None, state = 'inactive')
    Rule('CompI_formation', preCompI >> CompI(bDD=None, state = 'unmod'), Parameter('k2', 1))
    
    # --------------Complex I - RIP1 Modification-----------------
    # Complex I recruits RIP1
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='unmod'), 'bDD',[1e-6, 1e-3])
    #bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='ub'), 'bDD',[1e-6, 1e-3])
    
    # Complex I undergos modifications (e.g. Ubiquitination) and elicits NFkB signaling
    Rule('CompI_Ub', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'unmod')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), Parameter('k3', 1))
    #Rule('CompI_Ub2', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), Parameter('k4', 1))
    Rule('NFkB_signaling', CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub') + NFkB(bf=None), Parameter('k5', 1))
    
    # De-ubiquitination and degradation or release of RIP1 from complex I
    Rule('CompI_deUb', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub')>>CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,state='deub'), Parameter('k6',1e-3))
    
    Rule('RIP1_deg', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub') >> CompI(bDD=None, state='mod'), Parameter('k7',1))
    
    bind(CompI(bDD=None, state='mod'), 'bDD', RIP1(bDD=None, bRHIM = None,  state = 'deub'), 'bDD', [1e-10,1])
    
    # Receptor recycle
    Rule('TNFR1_recycle', CompI(bDD=None, state='mod') >> TRADD(bDD1=None, bDD2 = None, state='active') + TNFR1(blig = None, bDD = None, state =  'norm'), Parameter('k8',1e-3))
    
    # --------------RIP1 Ubiquitination---------------------------
    Rule('RIP1_Ub', RIP1(bDD=None, bRHIM = None, state='unmod')>> RIP1(bDD=None, bRHIM = None, state='ub'), Parameter('k9', 1e-7))

def SecondaryComplex_to_Bid_monomers():
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
    Parameter('BidK_0'  , 5.0e3) # molecules per cell
    
    alias_model_components()
    Initial(RIP3(bRHIM = None, state = 'unmod'), RIP3_0)   # RIP3
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
    Parameter('Ka_RIP1_FADD',   1e-7) # Biochemica et Biophysica Acta 1834(2013) 292-300
    Parameter('Kd_RIP1_FADD',   1e-8) # Biochemica et Biophysica Acta 1834(2013) 292-300
    alias_model_components()
    
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', TRADD(bDD1=None, state = 'active'), 'bDD1', [1e-6, 1e-3])
    bind(FADD(bDD = None), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'unmod'), 'bDD', [Ka_RIP1_FADD, Kd_RIP1_FADD])
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bDD', [1e-6, 1e-3])
    # For simplicity, I am neglecting the binary intereaction that occurs between proC8 and RIP1.
    # Binding of proC8 and c-flip to FADD is accomplished in CD95_to_Secondary complex. 

    #--------------RIP1 Truncation reactions-------------
    #---Truncation by C8---------------------------------
    RIP_CIIA_proC8 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% TRADD(bDD2 = None, bDD1 = ANY, state = 'active') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    RIP_CIIB_proC8 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    CIIA = TRADD(bDD2 = None, bDD1 = ANY, state = 'active') %  FADD(bDD=ANY, bDED1=None, bDED2=None)
    
    Rule('RIP1_truncation_CIIA', RIP_CIIA_proC8 >> CIIA + C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k11',1e-1))
    Rule('RIP1_truncation_CIIB', RIP_CIIB_proC8 >> FADD(bDD=None, bDED1=None, bDED2=None)+ C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k12', 1e-1))
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP1(bDD=None), 'bRHIM', 'state', 'unmod', 'trunc', [1e-6, 1e-3, 1e-1])

    #---Truncation by proC8:cFlip_L---------------------
    Riptosome_FADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%FADD(bDD=1, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Riptosome_TRADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%TRADD(bDD1=ANY, bDD2=1)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)

    Rule('RIP1_truncation_FADD', Riptosome_FADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k13', 1e-1))
    Rule('RIP1_truncation_TRADD', Riptosome_TRADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k14', 1e-1))
    
    # -------------RIP3 Binding Interactions----------------
    Ripto1_Flip_S = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Ripto2_Flip_S = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Necrosome1 = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=6, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 6, state = 'unmod')
    Necrosome2 = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=5, state='unmod') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 5, state = 'unmod')

    Rule('RIP3_binding1', Ripto1_Flip_S + RIP3(bRHIM= None, state = 'unmod') <> Necrosome1, Parameter('k15', 1e-6), Parameter('k16', 1e-3))
    Rule('RIP3_binding2', Ripto2_Flip_S + RIP3(bRHIM= None, state = 'unmod') <> Necrosome2, Parameter('k17', 1e-6), Parameter('k18', 1e-3))
    
    #RIP3 Truncation
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP3(), 'bRHIM', 'state', 'unmod', 'trunc', [1e-6, 1e-3, 1e-1])

    #-------------Bid Interactions--------------------------
    # Bid Phosphorylation and Truncation
    catalyze_state(BidK(), 'bf', Bid(), 'bf', 'state', 'U', 'po4', [1e-6, 1e-3, 1e-1])
    catalyze_state(C8(bf = None, state = 'A'), 'bf', Bid(), 'bf', 'state', 'U', 'T', [1.04e-5, 0.005, 0.1])

    # Bid-PO4 competing with RIP1 for binding to Complex II
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', Bid(bf = None, state = 'po4'), 'bf', [1e-6, 1e-3])
    # Bid-PO4 sequestering RIP1
    bind(RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bRHIM', Bid(bf = None, state = 'po4'), 'bf', [1e-6, 1e-3])

def SecondaryComplex_to_Bid_Alternate():
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
    Parameter('BidK_0'  , 5.0e3) # molecules per cell
    
    alias_model_components()
    Initial(RIP3(bRHIM = None, state = 'unmod'), RIP3_0)   # RIP3
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
    Parameter('Ka_RIP1_FADD',   1e-7) # Biochemica et Biophysica Acta 1834(2013) 292-300
    Parameter('Kd_RIP1_FADD',   1e-8) # Biochemica et Biophysica Acta 1834(2013) 292-300
    alias_model_components()
    
    #Assembling TRADD dependent Complex II
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', TRADD(bDD1=None, state = 'active'), 'bDD1', [1e-6, 1e-3])
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', RIP1(bDD = None, state = 'deub'), 'bDD', [1e-8, 1e-1])
    
    #Recruiting RIP1 to secondary complex and TRADD dependent Complex II
    bind(FADD(bDD = None, bDED1 = ANY, bDED2 =  ANY), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'unmod'), 'bDD', [Ka_RIP1_FADD, Kd_RIP1_FADD])
    bind(FADD(bDD = None, bDED1 = ANY, bDED2 =  ANY), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'deub'), 'bDD', [Ka_RIP1_FADD, Kd_RIP1_FADD])
    
    #bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bDD', [1e-6, 1e-1])
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'deub'), 'bDD', [1e-6, 1e-1])
    # For simplicity, I am neglecting the binary intereaction that occurs between proC8 and RIP1.
    # Binding of proC8 and c-flip to FADD is accomplished in CD95_to_Secondary complex.
    
    #--------------RIP1 Truncation reactions-------------
    #---Truncation by C8---------------------------------
    RIP_CIIA_proC8 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% TRADD(bDD2 = None, bDD1 = ANY, state = 'active') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    RIP_CIIA_proC8_alt = RIP1(bDD=ANY, bRHIM = None, state = 'deub')% TRADD(bDD2 = None, bDD1 = ANY, state = 'active') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    
    RIP_CIIB_proC8 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    RIP_CIIB_proC8_alt = RIP1(bDD=ANY, bRHIM = None, state = 'deub')% FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    
    CIIA = TRADD(bDD2 = None, bDD1 = ANY, state = 'active') %  FADD(bDD=ANY, bDED1=None, bDED2=None)
    
    Rule('RIP1_truncation_CIIA', RIP_CIIA_proC8 >> CIIA + C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k11',1e-1))
    Rule('RIP1_truncation_CIIA_alt', RIP_CIIA_proC8_alt >> CIIA + C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k11a',1e-6))
    
    Rule('RIP1_truncation_CIIB', RIP_CIIB_proC8 >> FADD(bDD=None, bDED1=None, bDED2=None)+ C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k12', 1e-1))
    Rule('RIP1_truncation_CIIB_alt', RIP_CIIB_proC8_alt >> FADD(bDD=None, bDED1=None, bDED2=None)+ C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k12a', 1e-6))
    
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP1(bDD=None), 'bRHIM', 'state', 'unmod', 'trunc', [1e-6, 1e-3, 1e-1])
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP1(bDD=None), 'bRHIM', 'state', 'deub', 'trunc', [1e-6, 1e-3, 1e-1])
    
    #---Truncation by proC8:cFlip_L---------------------
    Riptosome_FADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%FADD(bDD=1, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Riptosome_FADD_alt = RIP1(bDD=1, bRHIM = None, state = 'deub')%FADD(bDD=1, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    
    Riptosome_TRADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%TRADD(bDD1=ANY, bDD2=1)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Riptosome_TRADD_alt = RIP1(bDD=1, bRHIM = None, state = 'deub')%TRADD(bDD1=ANY, bDD2=1)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    
    Rule('RIP1_truncation_FADD', Riptosome_FADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k13', 1e-1))
    Rule('RIP1_truncation_FADD_alt', Riptosome_FADD_alt >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k13a', 1e-1))
    Rule('RIP1_truncation_TRADD', Riptosome_TRADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k14', 10))
    Rule('RIP1_truncation_TRADD_alt', Riptosome_TRADD_alt >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k14a', 10))
    
    # -------------RIP3 Binding Interactions----------------
    Ripto1_Flip_S = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Ripto2_Flip_S = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='unmod') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Necrosome1 = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=6, state='unmod') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 6, state = 'unmod')
    Necrosome2 = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=5, state='unmod') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 5, state = 'unmod')
    
    Ripto1_Flip_S_alt = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='deub') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Ripto2_Flip_S_alt = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=None, state='deub') % flip_S(bDED=ANY) % proC8(bDED=ANY)
    Necrosome1_alt = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=6, state='deub') % TRADD(bDD1=ANY, bDD2=ANY, state='active') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 6, state = 'unmod')
    Necrosome2_alt = FADD(bDD=ANY, bDED1=ANY, bDED2=ANY) % RIP1(bDD=ANY, bRHIM=5, state='deub') % flip_S(bDED=ANY) % proC8(bDED=ANY) % RIP3(bRHIM= 5, state = 'unmod')
    
    Rule('RIP3_binding1', Ripto1_Flip_S + RIP3(bRHIM= None, state = 'unmod') <> Necrosome1, Parameter('k15', 1e-6), Parameter('k16', 1e-3))
    Rule('RIP3_binding2', Ripto2_Flip_S + RIP3(bRHIM= None, state = 'unmod') <> Necrosome2, Parameter('k17', 1e-6), Parameter('k18', 1e-3))
    Rule('RIP3_binding1_alt', Ripto1_Flip_S_alt + RIP3(bRHIM= None, state = 'unmod') <> Necrosome1_alt, Parameter('k15a', 1e-6), Parameter('k16a', 1e-3))
    Rule('RIP3_binding2_alt', Ripto2_Flip_S_alt + RIP3(bRHIM= None, state = 'unmod') <> Necrosome2_alt, Parameter('k17a', 1e-6), Parameter('k18a', 1e-3))
    
    #RIP3 Truncation
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP3(), 'bRHIM', 'state', 'unmod', 'trunc', [1e-6, 1e-3, 1e-1])
    
    #-------------Bid Interactions--------------------------
    # Bid Phosphorylation and Truncation
    catalyze_state(BidK(), 'bf', Bid(), 'bf', 'state', 'U', 'po4', [1e-6, 1e-3, 1e-1])
    catalyze_state(C8(bf = None, state = 'A'), 'bf', Bid(), 'bf', 'state', 'U', 'T', [1.04e-5, 0.005, 0.1])
    
    # Bid-PO4 sequestering RIP1
    bind(RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bRHIM', Bid(bf = None, state = 'po4'), 'bf', [1e-6, 1e-3])
    bind(RIP1(bDD = None, bRHIM = None, state = 'deub'), 'bRHIM', Bid(bf = None, state = 'po4'), 'bf', [1e-6, 1e-3])

def rip1_to_parp():
    """Defines a connection between RIP1 and RIP3 and PARP activity observed in
    necroptosis. S. Joaun-Lanhouet et al. 2012 Obsered a requirement for RIP1
    RIP3 for PARP-1 activation. PARP-1 initiated TRAIL induced necroptosis in
    acidic extracellular conditions.

    Uses RIP1, RIP3 and MLKL monomers and their associated MLKL activation.
    """
    Monomer('MLKL', ['bRHIM', 'state'], {'state':['unmod', 'active', 'inactive']})
    Parameter('MLKL_0'  , 1.0e6) # molecules per cell
    alias_model_components()
    Initial(MLKL(bRHIM = None, state = 'unmod'), MLKL_0)   # MLKL
    
    Rule('Rip_PO4lation', RIP1(bRHIM=ANY, state = 'unmod')%RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'), Parameter('k19', 1e-1))
    Rule('Rip_PO4lation_alt', RIP1(bRHIM=ANY, state = 'deub')%RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'), Parameter('k19a', 1e-1))
    
    catalyze_state(RIP1(state='po4'), 'bPARP', MLKL(), 'bRHIM', 'state', 'unmod', 'active', [1e-6,1e-3, 1e-1])
    catalyze_state(MLKL(state='active'), 'bRHIM', MLKL(), 'bRHIM', 'state', 'unmod', 'active', [1e-7, 0.2, 0.01])
    
def C3_inhibits_MLKL():
    catalyze_state(C3(state='A'), 'bf', MLKL(), 'bRHIM', 'state','unmod','inactive', [1e-6, 1e-2, 1e-1])

def observables():
    Observable('Obs_TNFa', TNFa(blig =  None))
    Observable('Obs_Fas', Fas(blig = None))
    Observable('SecondaryComplex', FADD(bDD=None, bDED1 = ANY, bDED2 = ANY))
    Observable('Complex_IIA', TRADD(bDD1=ANY, bDD2=None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Riptosome1', RIP1(bDD = ANY, bRHIM = None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Riptosome2', RIP1(bDD = ANY, bRHIM = None)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Obs_RIP1_Ub', RIP1(state = 'ub'))
    Observable('RIP1_Bid', RIP1()%Bid())
    Observable('Bid_Riptosome1', Bid(bf= ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Bid_Riptosome2', Bid(bf= ANY)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Obs_cFlipS', flip_S())
    Observable('Obs_cFlipL', flip_L())
    Observable('Obs_NFkB', NFkB())
    Observable('Bid_Trunc', Bid(state='T'))
    Observable('Bid_PO4', Bid(state='po4'))
    Observable('RIP1_Trunc', RIP1(state='trunc'))
    Observable('RIP3_Trunc', RIP3(state='trunc'))
    Observable('Necrosome', RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'))
    Observable('Obs_proC8', proC8())
    Observable('Obs_C8', C8())
    Observable('Obs_C3ub', C3(state = 'ub'))
    Observable('Obs_C3', C3(state = 'A'))
    Observable('Obs_pC3', C3(state = 'pro'))
    #Observable('RIP1_nucl', RIP1(bDD = None, bRHIM=ANY, state = 'N')%RIP3(bRHIM=ANY, state = 'N'))
    Observable('Obs_cPARP', PARP(state='C'))
    Observable('Obs_PARP', PARP(state='U'))
    Observable('Obs_MLKL', MLKL(state = 'active'))