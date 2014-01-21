"""
    Overview
    ========
    
    PySB implementations of the apoptosis-necroptosis reaction model version 1.0
    (ANRM 1.0) originally published in [Irvin,NECRO2013]_.
    
    This model provides information about the dynamic biomolecular events that
    commit a cell to apoptosis or necroptosis in response to TNFa signalling.
    PARP1 cleavage or activity serves, in this case, as a marker for apoptosis 
    or necrosis respectively. An alternate marker for necroptosis
    
    This file contains functions that implement the parts of the extrinsic apoptosis
    pathway and activation of a necroptosis reporter proteins, PARP and MLKL. 
    
    The modules are organized into three overall sections:
    - ...

    These sections are further devided into...
    For receptor signalling and Bid activation:
    --...
    
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

# SECTION ONE: Fas Signalling
# ===========================
# This contains FasL binding and assembly of the DISC

def TNFa_to_ComplexI_Monomers():
    
    """ Declares TNFa, TNFR1, TRADD, RIP1, TRAF2, IAP, NKS (NFkB signaling complex),
    NFkB, CYLD and FADD. Binding sites are named after (when known) subdomains that
    facilitate binding. Otherwise the binding site will reflect what it binds with.
        
    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """
    Monomer('TNFa'  ,   ['brec'])           #TNFa
    Monomer('TNFR1' ,   ['blig', 'bDD'])    #TNFR1
    Monomer('TRADD' ,   ['bDD1', 'bDD2'])   #TRADD
    Monomer('RIP1'  ,   ['bDD', 'btraf', 'state'], {'state':['unmod', 'ub']})   #RIP1
    Monomer('TRAF'  ,   ['brip', 'bciap', 'zfinger'])   #TRAF2
    Monomer('cIAP'  ,   ['btraf'])          #cIAP 1 and 2
    Monomer('NSC'   ,   ['bnfkb'])          #Recruitement of LUBAC commits complex I to the NFkB activation pathway. This pathway is approximated by a single activation step which is carried out by a NFkB signaling complex (NSC). 
    Monomer('NFkB'  ,   ['bnsc','state'], {'state':['I', 'A']})       #NFkB
    Monomer('CYLD'  ,   ['btraf'])  #CYLD
    Monomer('FADD', ['bDD', 'bDED1','bDED2'])    #FADD
    alias_model_components()

def TNFa_to_ComplexI_Initials():
    Parameter('TNFa_0'  ,  6000) # 6000 corresponds to 100ng/ml TNFa
    Parameter('TNFR1_0' ,  4800) # 4800 receptors per cell
    Parameter('TRADD_0' ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('RIP1_0'  , 12044) # molecules per cell (arbitrarily assigned) 12044
    Parameter('TRAF_0'  ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('cIAP_0'  ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('NSC_0'   ,     0) # complexes per cell
    Parameter('NFkB_0'    , 50000) # molecules per cell (arbitrarily assigned) 50000
    Parameter('CYLD_0'  ,  9000) # molecules per cell
    Parameter('FADD_0'  ,   8030) # 8300 molecules per cell (J Immunol 2005)
    alias_model_components()
    
    Initial(TNFa(brec=None), TNFa_0)       
    Initial(TNFR1(blig=None, bDD=None), TNFR1_0)
    Initial(TRADD(bDD1=None, bDD2=None), TRADD_0)
    Initial(RIP1(bDD=None, btraf=None, state = 'unmod'), RIP1_0)
    Initial(TRAF(brip=None, bciap=None, zfinger=None), TRAF_0)
    Initial(cIAP(btraf=None), cIAP_0)
    Initial(NSC(bnfkb=None), NSC_0)
    Initial(NFkB(bnsc=None, state = 'I'), NFkB_0)
    Initial(CYLD(btraf=None),CYLD_0)
    Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)

def TNFa_to_ComplexI():
    """Reaction network that produces Complex I and activates NFkB. Instead of including conversion step. 
    """
    bind(TNFa(brec = None),'brec', TNFR1(blig = None), 'blig', [1e-6, 1e-3])
    bind(TNFR1(blig = ANY, bDD = None), 'bDD', TRADD(bDD1 = None, bDD2 = None), 'bDD1', [1e-6, 1e-3])
    bind(TNFR1(blig = ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = None), 'bDD2', RIP1(bDD = None, btraf = None, state =  'unmod'), 'bDD', [1e-6, 1e-3])
    bind(TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY, btraf = None, state =  'unmod'), 'btraf', TRAF(brip=None), 'brip', [1e-6, 1e-3])
    bind(TRAF(bciap =  None), 'bciap', cIAP(btraf = None), 'btraf', [1e-6, 1e-3])
    bind(TRAF(zfinger =  None), 'zfinger', CYLD(btraf = None), 'btraf', [1e-6, 1e-3])
    
    #Make these aliases.
    ComplexI = TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)
    
    rule('RIP1_Ubiquitination', ComplexI %RIP1(bDD = ANY, btraf = None, state =  'unmod')%TRAF(brip = ANY, bciap = ANY) >> ComplexI %RIP1(bDD = ANY, btraf = None, state =  'ub')%TRAF(brip = ANY, bciap = ANY), Parameter('RIP1_ubiq_kc', 1e-1))
    rule('RIP1_Deubiquitination', ComplexI%RIP1(bDD = ANY, btraf = None, state =  'ub')%TRAF(brip = ANY, bciap = ANY) >> ComplexI%RIP1(bDD = ANY, btraf = None, state =  'unmod')%TRAF(brip = ANY, zfinger = ANY), Parameter('RIP1_deub_kc', 1e-1))
    rule('Establish_NFkB_Signaling_Complex', ComplexI%RIP1(bDD = ANY, btraf = None, state =  'ub')%TRAF(brip = ANY, bciap = ANY) >> NSC(bnfkb=None), Parameter('NSC_esta_kc', 1e-1))
    rule('RIP1_Degradation', ComplexI%RIP1(bDD = ANY, btraf = None, state =  'ub')%TRAF(brip = ANY) >> ComplexI + TRAF(brip = None), Parameter('RIP1_degr_kc', 1e-1))

def CompI_TRADD_RIP1_Dissociation():
    """Some believe that TRADD, RIP1 and others are released into the cytoplasm after RIP1 deubiquitination or after TNFR endocytosis. This mechanism presents a problem: The RIP1 released into the cytoplasm is in theory the same as RIP1 that existed in the cytoplasm prior to TNF stimulation. Meaning apoptosis whould spontaneously occur. Dissociation of TRADD:RIP1 from complex I is a hypothetical pathway that could distinguish pre- and post- TNF RIP1. It could also serve as a mechanism for FADD independent necrosome formation. This hypothesis is supported by a suggested mechanism, to explain signal transduction from TNFa to Riptosome formation[1]. 
        
        1. Dickens, LS., IR Powley, MA Hughes, M MacFarlane, Exp. Cell. Res. 318 (2012) 1269-1277"""
    
    bind(TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = None), 'bDD', TRADD(bDD1 = None, bDD2 = ANY)%RIP1(bDD = ANY, state =  'unmod'), 'bDD1', [1e-6, 1e-3])

def CompII_Hypothesis_1():
    """Michaeu and Tschopp, showed that FADD is absent from Complex I[1]. The presence of TRADD in Complex II is debated [2]. As explained in CompI_TRADD_RIP1_Dissociation(), simply releasing RIP1 from Complex I so that it can bind FADD creates a condition where cell death spontaneously occurs. Here, I hypothesize that FADD transiently associates with Complex I in order to retrieve RIP1 from Complex I. This association, being transient, would have been difficult to detect.
        
        1. Olivier Micheau, Jurg Tschopp, Induction of TNF Receptor I-Mediated Apoptosis via Two Sequential Signaling Complexes, Cell, Volume 114, Issue 2, 25 July 2003, Pages 181-190
        
        2. Dickens, LS., IR Powley, MA Hughes, M MacFarlane, "The 'complexities' of life and death: Death receptor signalling platforms" Exp. Cell. Res. 318 (2012) 1269-1277"""

    TNFR1(bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY)%TRAF2(brip = ANY) + FADD(bDD = None)>> TNFR1(bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = None) + RIP1(bDD1 = ANY)%FADD(bDD = ANY) + TRAF2(brip = None)
    
    TNFR1(bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY)+ FADD(bDD = None)>> TNFR1(bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = None) + RIP1(bDD1 = ANY)%FADD(bDD = ANY)
    


def FADD_to_C8():
    """Defines the interactions of procaspase 8 and/or cFlip recruitment to FADD as per
    ANRM 1.0.
        
    Uses FADD, proC8 and flip_L and flip_S monomers and parameters to generate the rules.
    
    FADD is a scaffold protein that gains the ability to recruit cFlip and procaspase 8 when
    it's death domain DD is occupied. Fas, Rip1 and TRADD are capable of binding FADD. 
    """

    # ================================================
    # c-Flip and procaspase 8 recruitment to FADD
    # ------------------------------------------------
    #   X:FADD:proC8 + proC8  <-> X:FADD:proC8:proC8
    #   X:FADD:proC8 + flip_L <-> X:FADD:proC8:flip_L
    #   X:FADD:proC8 + flip_S <-> X:FADD:proC8:flip_S
    #
    #   X = RIP1, TRADD or Fas

    # ----------procaspase8 and cFlip recruitment-----
    bind(FADD(bDD = ANY, bDED1 = None, bDED2 = None),'bDED1', proC8(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    
    # procaspase 8 dimerization
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', proC8(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    
    # For simplicity allow proC8 to bind FADD before any c-Flip do.
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', flip_L(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    bind(FADD(bDD = ANY, bDED2 = None, bDED1 = ANY), 'bDED2', flip_S(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)

    # Caspase 8 activation
    Rule('C8_activation', FADD(bDED2 = ANY, bDED1 = ANY)%proC8(bDED = ANY)%proC8(bDED = ANY) >> FADD(bDED1 = None, bDED2 = None) + C8(bf = None, state = 'A'), Parameter('kc_c8', 1))

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
    """One paper mentioned that cFlip present in the secondary complex inhibits necrosis by allowing procaspase8 mediated cleavage of RIP1. But, in the presence of cFlip, procaspase 8 cannot adaquately cleave Bid to initiate apoptosis. Its was questioned by the fact that Bid has a higher affinity for caspase 8 than RIP1. However, Dr. Zinkel later found that Bid is absent from the Secondary complex. This hypothetical interaction between cFlip and Bid is presented here:"""
    Rule('Bid_to_SecComp', FADD(bDD = None, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED = ANY)%proC8(bDED = ANY) + Bid(bf=None, state='U') <> FADD(bDD = 1, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED = ANY)%proC8(bDED = ANY)%Bid(bf=1, state='U'), Parameter('Bid_Flip_Bindf', 1e-6),Parameter('Bid_Flip_Bindr', 1e-1))
    
    Rule('tBid_SlowFormation', FADD(bDD = ANY, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED = ANY)%proC8(bDED = ANY)%Bid(bf=ANY, state='U') >> FADD(bDD = None, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED = ANY)%proC8(bDED = ANY) + Bid(bf=None, state='T'), Parameter('Bid_Flip_cat', 1e-1))

def TNFR1_to_ComplexII_monomers():
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

def TNFR1_to_ComplexII_Initials():
    Parameter('TNFa_0'  ,  6000) # 6000 corresponds to 100ng/ml TNFa
    Parameter('TNFR1_0' ,  4800) # 4800 receptors per cell
    Parameter('TRADD_0' ,  9000) # molecules per cell (arbitrarily assigned) 9000
    Parameter('CompI_0' ,     0) # complexes per cell
    Parameter('NFkB_0'  ,     0) # molecules per cell
    alias_model_components()
    
    Initial(TNFa(blig=None), TNFa_0)                                 # TNFa Ligand
    Initial(TNFR1(blig=None, bDD=None, state='norm'), TNFR1_0)       # TNFR1
    Initial(TRADD(bDD1=None, bDD2=None, state='inactive'), TRADD_0)  # TRADD
    Initial(CompI(bDD=None, state='unmod'), CompI_0)                 # Complex I
    Initial(NFkB(bf=None), NFkB_0)

def TNFR1_to_ComplexI():
    """TNFR1_recruits TRADD and TRADD which serves as a scaffold for RIP1, ubiquitinating complexes and d-ubiquitinating complexes. RIP1 may also bind to TNFR1 in the absence of TRADD.
        http://www.nature.com/ni/journal/v9/n9/fig_tab/ni0908-1015_F1.html
        We are only considering the case where complex I formation requires TRADD.
        
        Following TRADD recruitment to TNFR1, ubiqutinating and de-ubiquitinating molecules bind modify protiens in Complex I (possibly including TRADD). The hypothesized modification of TRADD in complex I is included in this module."""
    
    # =========================================
    # TNFR1 ligation, formation of Complex I and release of RIP1 and TRADD rules
    # -----------------------------------------
    #   TNFa+ TNFR1 <-> TNFa:TNFR1
    #   TNFa:TNFR1 + TRADD <-> TNFa:TNFR1:TRADD >> CompI
    # ------------------------------------------
    
    # -------------Complex I assembly----------------
    bind(TNFa(blig=None), 'blig', TNFR1(blig=None, bDD=None, state='norm'), 'blig', [2e-7, 1e-4])
    bind(TNFR1(blig = ANY, bDD = None, state =  'norm'), 'bDD', TRADD(bDD1=None, bDD2=None, state='inactive'), 'bDD1', [1e-6, 0.2])
    preCompI = TNFa(blig=ANY)%TNFR1(blig=ANY, bDD=ANY, state = 'norm')%TRADD(bDD1 = ANY, bDD2=None, state = 'inactive')
    Rule('CompI_formation', preCompI >> CompI(bDD=None, state = 'unmod'), Parameter('k2', 1))

def RIP1_Hypothesis_1():
    """Deubiquitinated RIP1, released from Complex I, can bind FADD and initiate apoptosis and necrosis (the literature). Yet, there is nothing explicitely stated in the literature to distinguish deubiquitinated from unmodifed RIP1. If RIP1-unmod = RIP1-deub (i.e. RIP1 Hypothesis 1), then unstimulated cell will undergo apoptosis. There needs to be a starting point that is not RIP1-unmod. Papers show cIAP expression inhibits apoptosis by ubiquitinating RIP1. If we let the initial state of the cell be primarily or ALL RIP1-ub, then we can have a cell that reponses to pro-apoptotic stimuli as demonstrated in literature."""
    Parameter('RIP1_0'  , 12044) # molecules per cell (arbitrarily assigned) 12044
    alias_model_components()
    Initial(RIP1(bDD=None, bRHIM = None, bPARP = None, state = 'ub'), RIP1_0)   # RIP1
    
    # =============================
    # Complex I - RIP1 Modification
    # -----------------------------
    #   CompI + RIP1-ub <> CompI:RIP1-ub >> [active]CompI:RIP1-Ub
    #   [active]CompI:RIP1-Ub >> [active]CompI:RIP1
    #   [active]CompI:RIP1-Ub <> [active]CompI + RIP1
    #   [active]CompI:RIP1-Ub >> [active]CompI
    # -----------------------------
    
    # Complex I recruits RIP1
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='ub'), 'bDD',[1e-6, 1e-3])
    Rule('CompI_Ub_2', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), Parameter('k3', 1))

    # De-ubiquitination and degradation or release of RIP1 from complex I
    Rule('CompI_deUb', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub')>>CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,state='unmod'), Parameter('k6',1e-3))
    bind(CompI(bDD=None, state='mod'), 'bDD', RIP1(bDD=None, bRHIM = None,  state = 'unmod'), 'bDD', [1e-10,1])
    Rule('RIP1_deg', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub') >> CompI(bDD=None, state='mod'), Parameter('k7',1))

def RIP1_Hypothesis_2():
    """Deubiquitinated RIP1, released from Complex I, can bind FADD and initiate apoptosis and necrosis (the literature). Yet, there is nothing explicitely stated in the literature to distinguish deubiquitinated from unmodifed RIP1. Uniprot reports a binding reaction between RIP1 and FADD. If RIP1-unmod =! RIP1-deub (i.e. RIP1 Hypothesis 2), then a RIP1-deub-FADD binding interaction could facilitate the apoptotic/necrotic signaling that occurs downstream of TNFR1 without creating a state that is """

    Parameter('RIP1_0'  , 12044) # molecules per cell (arbitrarily assigned) 12044
    alias_model_components()
    Initial(RIP1(bDD=None, bRHIM = None, bPARP = None, state = 'unmod'), RIP1_0)   # RIP1
    # =============================
    # Complex I - RIP1 Modification
    # -----------------------------
    #   CompI + RIP1 <> CompI:RIP1 >> [active]CompI:RIP1-Ub
    #   [active]CompI:RIP1-Ub >> [active]CompI:RIP1-deub
    #   [active]CompI:RIP1-deub <> [active]CompI + RIP1-deub
    #   [active]CompI:RIP1-deub >> [active]CompI
    # -----------------------------
    
    # Complex I recruits RIP1
    bind(CompI(bDD=None, state = 'unmod'), 'bDD', RIP1(bDD=None, bRHIM = None, state='unmod'), 'bDD',[1e-6, 1e-3])
    Rule('CompI_Ub_1', CompI(bDD=ANY, state = 'unmod')%RIP1(bDD=ANY,bRHIM=None, state = 'unmod')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub'), Parameter('k3', 1))
    
    # De-ubiquitination and degradation or release of RIP1 from complex I
    Rule('CompI_deUb', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='ub')>>CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,state='deub'), Parameter('k6',1e-3))
    
    Rule('RIP1_deg', CompI(bDD=ANY, state='mod')%RIP1(bDD=ANY, bRHIM=None,  state='deub') >> CompI(bDD=None, state='mod'), Parameter('k7',1))
    
    bind(CompI(bDD=None, state='mod'), 'bDD', RIP1(bDD=None, bRHIM = None,  state = 'deub'), 'bDD', [1e-10,1])

def ComplexI_to_NFkB():
    Rule('NFkB_signaling', CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub')>> CompI(bDD=ANY, state = 'mod')%RIP1(bDD=ANY,bRHIM=None, state = 'ub') + NFkB(bf=None), Parameter('k5', 1))
    
    # Receptor recycle
    Rule('TNFR1_recycle', CompI(bDD=None, state='mod') >> TRADD(bDD1=None, bDD2 = None, state='active') + TNFR1(blig = None, bDD = None, state =  'norm'), Parameter('k8',1e-3))

def SecondaryComplex_to_Bid_monomers():
    Monomer('BidK', ['bf']) #unknown Bid-kinase
    Monomer('RIP3', ['bRHIM', 'state'], {'state':['unmod', 'po4', 'trunc', 'N']})
    alias_model_components()

def SecondaryComplex_to_Bid_initials():
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

def ComplexII_Hypothesis_1():
    
    # ==============================================================
    # Initializataion of Complex II
    # --------------------------------------------------------------
    #   FADD + TRADD[active] <-> FADD:TRADD[active]
    #   FADD + RIP1 <-> FADD:RIP1
    #   TRADD + RIP1 <-> TRADD:RIP1
    # --------------------------------------------------------------

    # -------------Assembling Complex II-----------------
    Parameter('Ka_RIP1_FADD',   1e-7) # Biochemica et Biophysica Acta 1834(2013) 292-300
    Parameter('Kd_RIP1_FADD',   1e-8) # Biochemica et Biophysica Acta 1834(2013) 292-300
    alias_model_components()
    
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', TRADD(bDD1=None, state = 'active'), 'bDD1', [1e-6, 1e-3])
    bind(FADD(bDD = None), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'unmod'), 'bDD', [Ka_RIP1_FADD, Kd_RIP1_FADD])
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bDD', [1e-6, 1e-3])

def ComplexII_Hypothesis_2():
    """The additional RIP1 state provided in RIP1_Hypothesis_2() creates a problem for FasL signaling. One hypothesis we have is that the signaling networks responsible for apoptosis and necrosis are similar for Fas and TNF sensitive cells. Under this assumption, RIP1-deub might bind FADD. FADD would then recruit procaspase 8 and/or cFlip to form TNFR1 dependent complext II. At this point RIP1-deub may dissociate and be replaced by RIP1-unmod (FADD-cFlip-procaspase8 == Fas dependent secondary complex and is capable of binding RIP1-unmod). This effect did not seem to alter the dynamics much. TRADD's ability to recruit RIP-unmod, however did!"""
        
    #Assembling TRADD dependent Complex II
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', TRADD(bDD1=None, state = 'active'), 'bDD1', [1e-6, 1e-3])
    bind(FADD(bDD = None, bDED1 = None, bDED2 = None), 'bDD', RIP1(bDD = None, state = 'deub'), 'bDD', [1e-8, 1e-1])
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'deub'), 'bDD', [1e-6, 1e-1])

def TRADD_RIP1unmod_Hypothesis_2():
    bind(TRADD(bDD2 = None, state = 'active'),'bDD2', RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bDD', [1e-6, 1e-1])
    
def RIP1_to_SecondaryComplex():
    # -------------Assembling Complex II-----------------
    Parameter('Ka2_RIP1_FADD',   1e-7) # Biochemica et Biophysica Acta 1834(2013) 292-300
    Parameter('Kd2_RIP1_FADD',   1e-8) # Biochemica et Biophysica Acta 1834(2013) 292-300
    alias_model_components()
    # ==============================================================
    # Recruitment of Rip1 to the Secondary Complex
    # --------------------------------------------------------------
    #   FADD:x:x + RIP1unmod <-> FADD:x:x:RIP1
    #   FADD:x:x + RIP1deub  <-> FADD:x:x:RIP1  <<if Rip1 disociates from the complex II
    # --------------------------------------------------------------
    
    #Recruiting RIP1 to secondary complex and TRADD dependent Complex II
    bind(FADD(bDD = None, bDED1 = ANY, bDED2 =  ANY), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'unmod'), 'bDD', [Ka2_RIP1_FADD, Kd2_RIP1_FADD])
    bind(FADD(bDD = None, bDED1 = ANY, bDED2 =  ANY), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'deub'), 'bDD', [Ka2_RIP1_FADD, Kd2_RIP1_FADD])
    bind(FADD(bDD = None, bDED1 = ANY, bDED2 =  ANY), 'bDD', RIP1(bDD=None, bRHIM = None, state = 'ub'), 'bDD', [Ka2_RIP1_FADD, Kd2_RIP1_FADD])

def RIP1_deubiqutination_Hypothesis_1():
    """Necrosomes assemble around a deubiquitinated RIP1. In hypothesis1, Fas signaling does not deubiquitinate RIP1-ub. Without RIP1 deubiquitination the necrosome cannot form. This module imposes a RIP1-deubiquitinatoin step that would aloow necrosome formation in Fas sensitive cells"""
    Rule('RIP1ub_to_unmod', FADD(bDD = ANY, bDED1 = ANY, bDED2 =  ANY)% RIP1(bDD=ANY, bRHIM = None, state = 'ub') >> FADD(bDD = ANY, bDED1 = ANY, bDED2 =  ANY)% RIP1(bDD=ANY, bRHIM = None, state = 'unmod'), Parameter('kc_rip1', 1e-3))

def RIP1_truncation():
    """RIP1 Truncation as it occurs in the secondary complex and/or complex II. When RIP1 shares
        the complex with cFlip-S, RIP1 is stabilized and can bind RIP3 to form the necrosome."""
    
    # ==============================================================
    # Rip1 truncation by the Secondary Complex
    # --------------------------------------------------------------
    #   Rip1:TRADD:FADD:proC8:proC8 >> TRADD:FADD + Rip1-trunc + C8
    #   RIP1:FADD:proC8:proC8 >> FADD + Rip1-trunc + C8
    #
    #   Rip1:TRADD:FADD:proC8:cFlip >> TRADD:FADD:proC8:cFlip + Rip1-trunc
    #   RIP1:FADD:proC8:cFlip >> FADD:proC8:cFlip + Rip1-trunc
    #
    #   RIP1%ProC8%ProC8(in a complex) >> RIP1[trunc] + C8 + (remains of the complex)
    #   RIP1%ProC8%cFlip[L](in a complex) >> RIP1[trunc] + remains of the complex)
    #   RIP1%cFlip[S](in a complex) + RIP3 >> RIP1:RIP3(in a complex, i.e. necrosome)
    # --------------------------------------------------------------
    
    #--------------RIP1 Truncation reactions-------------
    #---Truncation by C8---------------------------------
    RIP_CIIA_proC8_1 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% TRADD(bDD2 = None, bDD1 = ANY, state = 'active') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    RIP_CIIA_proC8_2 = RIP1(bDD=ANY, bRHIM = None, state = 'deub')% TRADD(bDD2 = None, bDD1 = ANY, state = 'active') % FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    
    RIP_CIIB_proC8_1 = RIP1(bDD=ANY, bRHIM = None, state = 'unmod')% FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    RIP_CIIB_proC8_2 = RIP1(bDD=ANY, bRHIM = None, state = 'deub')% FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED=ANY)%proC8(bDED=ANY)
    CIIA = TRADD(bDD2 = None, bDD1 = ANY, state = 'active') %  FADD(bDD=ANY, bDED1=None, bDED2=None)
    
    Rule('RIP1_truncation_CIIA_1', RIP_CIIA_proC8_1 >> CIIA + C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k11',1e-1))
    Rule('RIP1_truncation_CIIA_2', RIP_CIIA_proC8_2 >> CIIA + C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k11a',1e-6))
    
    Rule('RIP1_truncation_CIIB_1', RIP_CIIB_proC8_1 >> FADD(bDD=None, bDED1=None, bDED2=None)+ C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k12', 1e-1))
    Rule('RIP1_truncation_CIIB_2', RIP_CIIB_proC8_2 >> FADD(bDD=None, bDED1=None, bDED2=None)+ C8(bf = None, state = 'A') + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k12a', 1e-6))
    
    #---Truncation by proC8:cFlip_L---------------------
    Riptosome_FADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%FADD(bDD=1, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Riptosome_FADD_alt = RIP1(bDD=1, bRHIM = None, state = 'deub')%FADD(bDD=1, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    
    Riptosome_TRADD = RIP1(bDD=1, bRHIM = None, state = 'unmod')%TRADD(bDD1=ANY, bDD2=1)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    Riptosome_TRADD_alt = RIP1(bDD=1, bRHIM = None, state = 'deub')%TRADD(bDD1=ANY, bDD2=1)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY)
    
    Rule('RIP1_truncation_FADD_1', Riptosome_FADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k13', 1e-1))
    Rule('RIP1_truncation_FADD_2', Riptosome_FADD_alt >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k13a', 1e-1))
    Rule('RIP1_truncation_TRADD_1', Riptosome_TRADD >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k14', 10))
    Rule('RIP1_truncation_TRADD_2', Riptosome_TRADD_alt >> FADD(bDD=None, bDED1=ANY, bDED2=ANY)%proC8(bDED = ANY)%flip_L(bDED = ANY) + RIP1(bDD=None, bRHIM = None, state = 'trunc'), Parameter('k14a', 10))
    
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
    
def Bid_Hypothesis():
    """The Zinkel lab found evidence of Bid mediated inhibition of necroptosis. Dr. Zinkel proposed that Bid has a third state (e.g. phosporylated) that can inhibit the necrosome by sequestering RIP3. We presume that Bid-po4 cannot be truncated. There are 3 phosphorylation sites on Bid and evidence that BId-po4 is resistant to truncation exists"""
    
    #-------------Bid Interactions--------------------------
    # Bid Phosphorylation and Truncation
    catalyze_state(BidK(), 'bf', Bid(), 'bf', 'state', 'U', 'po4', [1e-6, 1e-3, 1e-1])
    
    # Bid-PO4 sequestering RIP1
    bind(RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bRHIM', Bid(bf = None, state = 'po4'), 'bf', [1e-6, 1e-3])
    bind(RIP1(bDD = None, bRHIM = None, state = 'deub'), 'bRHIM', Bid(bf = None, state = 'po4'), 'bf', [1e-6, 1e-3])

def Bid_RIP_recruits_proC8():
    Rule('Bid_recruits_1', RIP1(bDD = None, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') + proC8(bDED = None) <> RIP1(bDD = 1, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') % proC8(bDED = 1), Parameter('kbid1', 1e-6), Parameter('kbid2', 1e-3))

    Rule('Bid_recruits_2', RIP1(bDD = None, bRHIM =  ANY, state = 'deub') % Bid(bf = ANY, state = 'po4') + proC8(bDED = None) <> RIP1(bDD = 1, bRHIM =  ANY, state = 'deub') % Bid(bf = ANY, state = 'po4') % proC8(bDED = 1), Parameter('kbid3', 1e-6), Parameter('kbid4', 1e-3))

def Bid_RIP_proC8_truncation():
    Rule('Bid_trunc_RIP_1', RIP1(bDD = ANY, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') % proC8(bDED = ANY) >> Bid(bf = None, state = 'po4') + RIP1(bDD = None, bRHIM = None, state = 'trunc') + proC8(bDED = None), Parameter('kbid5', 1e-1))
    
    Rule('Bid_trunc_RIP_2', RIP1(bDD = ANY, bRHIM =  ANY, state = 'deub') % Bid(bf = ANY, state = 'po4') % proC8(bDED = ANY) >> Bid(bf = None, state = 'po4') + RIP1(bDD = None, bRHIM = None, state = 'trunc')+ proC8(bDED = None), Parameter('kbid6', 1e-1))

def Bid_RIP_proC8_to_NFkB():
    Rule('Bid_RIP_NFkB', RIP1(bDD = ANY, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') % proC8(bDED = ANY) >>RIP1(bDD = ANY, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') % proC8(bDED = ANY) + NFkB(bf=None), Parameter('kbid7', 1e-1))

def Bidpo4_to_tBid_Hypothesis():
    """Here we allow truncation of Bid_po4"""
    catalyze_state(C8(bf = None, state = 'A'), 'bf', Bid(), 'bf', 'state', 'po4', 'T', [1.04e-5, 0.005, 0.1])

def C8_catalyzed_truncations():
    
    #   RIP1 + C8 <-> RIP1:C8 >> RIP1[trunc] + C8
    #   RIP3 + C8 <-> RIP3:C8 >> RIP3[trunc] + C8
    #   Bid + C8 <-> Bid:C8 >> Bid[trunc] + C8
    
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP1(bDD=None), 'bRHIM', 'state', 'unmod', 'trunc', [1e-6, 1e-3, 1e-1])
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP1(bDD=None), 'bRHIM', 'state', 'deub', 'trunc', [1e-6, 1e-3, 1e-1])
    
    #RIP3 Truncation
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP3(), 'bRHIM', 'state', 'unmod', 'trunc', [1e-6, 1e-3, 1e-1])

    #Bid Truncation
    catalyze_state(C8(bf = None, state = 'A'), 'bf', Bid(), 'bf', 'state', 'U', 'T', [1.04e-5, 0.005, 0.1])

def rip1_to_MLKL_monmers():
    Monomer('MLKL', ['bRHIM', 'state'], {'state':['unmod', 'active', 'inactive']})
    
def rip1_to_MLKL_initials():
    Parameter('MLKL_0'  , 1.0e6) # molecules per cell
    alias_model_components()
    Initial(MLKL(bRHIM = None, state = 'unmod'), MLKL_0)   # MLKL
    
def rip1_to_MLKL():
    """MLKL has recently been established as a substrate of the necrosome 3. And is considered an effector of necroptosis. However, it has not been established whether MLKL is de-activated in apoptosis. (Expected downstream proteins (i.e. Drp) are active in apoptosis as well as necrosis).
    
        Uses RIP1, RIP3 and MLKL monomers and their associated MLKL activation.
        """
    Rule('Rip_PO4lation', RIP1(bRHIM=ANY, state = 'unmod')%RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'), Parameter('k19', 1e-1))
    Rule('Rip_PO4lation_alt', RIP1(bRHIM=ANY, state = 'deub')%RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'), Parameter('k19a', 1e-1))
    
    catalyze_state(RIP1(state='po4'), 'bPARP', MLKL(), 'bRHIM', 'state', 'unmod', 'active', [1e-6,1e-3, 1e-1])
    catalyze_state(MLKL(state='active'), 'bRHIM', MLKL(), 'bRHIM', 'state', 'unmod', 'active', [1e-7, 0.2, 0.01])

def C3_inhibits_MLKL():
    """it has not been established whether MLKL is de-activated in apoptosis. But inorder to make a cellular decision you have to have one choice, once selected, inhibit all of the alternative pathways. This is a hypothetical apoptosis mediated inhibition of necrosis effector, MLKL."""
    catalyze_state(C3(state='A'), 'bf', MLKL(), 'bRHIM', 'state','unmod','inactive', [1e-6, 1e-2, 1e1])

def Momomers_zVad_to_C8():
    Monomer('zVad', ['bC8'])       # apoptosis inhibitor
    alias_model_components()

def Initials_zVad_to_C8():
    Parameter('zVad_0'   ,  0) # 20uM zVad converts to 9.6e6  molecules per cell if the cell volume = 8e-13L.
    alias_model_components()
    
    Initial(zVad(bC8=None), zVad_0)

def zVad_to_C8():
    bind(zVad(bC8 = None), 'bC8', C8(bf = None, state = 'A'), 'bf', [1e-6, 1e-3])
    bind(zVad(bC8 = None), 'bC8', C3(bf = None, state = 'A'), 'bf', [1e-6, 1e-3])
    
    Rule('zVad_C8', zVad(bC8 = ANY)%C8(bf = ANY, state = 'A') >> zVad(bC8 = ANY)%C8(bf = ANY, state = 'I'), Parameter('kzVad1', 1e-1))
    Rule('zVad_C3', zVad(bC8 = ANY)%C3(bf = ANY, state = 'A') >> zVad(bC8 = ANY)%C3(bf = ANY, state = 'I'),Parameter('kzVad2', 1e-1))
         
def observables():
    Observable('Obs_TNFa', TNFa(blig =  None))
    Observable('Obs_Fas', Fas(blig = None))
    Observable('SecondaryComplex', FADD(bDD=None, bDED1 = ANY, bDED2 = ANY))
    Observable('Complex_IIA', TRADD(bDD1=ANY, bDD2=None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Riptosome1', RIP1(bDD = ANY, bRHIM = None)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Riptosome2', RIP1(bDD = ANY, bRHIM = None)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Obs_RIP1_Ub', RIP1(state = 'ub'))
    Observable('RIP1_Bid', RIP1()%Bid())
    Observable('secondary_Complex', FADD(bDD=None, bDED1=ANY, bDED2=ANY))
    Observable('Bid_Riptosome1', Bid(bf= ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Bid_Riptosome2', Bid(bf= ANY)%TRADD(bDD1=ANY, bDD2=ANY)%FADD(bDD=ANY, bDED1=ANY, bDED2=ANY))
    Observable('Obs_cFlipS', flip_S())
    Observable('Obs_cFlipL', flip_L())
    Observable('Obs_NFkB', NFkB())
    Observable('Bid_Trunc', Bid(state='T'))
    Observable('Bid_PO4', Bid(state='po4'))
    Observable('Obs_RIP1', RIP1(bDD = None))
    Observable('RIP1_Trunc', RIP1(state='trunc'))
    Observable('RIP3_Trunc', RIP3(state='trunc'))
    Observable('Necrosome', RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'))
    Observable('Obs_proC8', proC8())
    Observable('Obs_C8', C8())
    Observable('Obs_C3ub', C3(state = 'ub'))
    Observable('Obs_C3', C3(state = 'A'))
    Observable('Obs_pC3', C3(state = 'pro'))
    Observable('RIP1_FADD', FADD(bDD = ANY)%RIP1(bDD=ANY))
    #Observable('RIP1_nucl', RIP1(bDD = None, bRHIM=ANY, state = 'N')%RIP3(bRHIM=ANY, state = 'N'))
    Observable('Obs_cPARP', PARP(state='C'))
    Observable('Obs_PARP', PARP(state='U'))
    Observable('Obs_MLKL', MLKL(state = 'active'))
    Observable('Obs_CytoC', CytoC(state='A'))
