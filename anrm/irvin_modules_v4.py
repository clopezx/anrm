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
    Monomer('RIP1'  ,   ['bDD', 'btraf', 'bRHIM', 'bMLKL', 'state'], {'state':['unmod', 'ub', 'po4', 'trunc']})   #RIP1
    Monomer('TRAF'  ,   ['brip', 'bciap', 'zfinger'])   #TRAF2
    Monomer('cIAP'  ,   ['btraf'])          #cIAP 1 and 2
    Monomer('NSC'   ,   ['bnfkb'])          #Recruitement of LUBAC commits complex I to the NFkB activation pathway. This pathway is approximated by a single activation step which is carried out by a NFkB signaling complex (NSC). 
    Monomer('NFkB'  ,   ['bnsc','state'], {'state':['I', 'A']})       #NFkB
    Monomer('CYLD'  ,   ['btraf', 'state'], {'state':['U','T']})  #CYLD
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
    Initial(RIP1(bDD=None, btraf=None, bRHIM = None, bMLKL = None, state = 'unmod'), RIP1_0)
    Initial(TRAF(brip=None, bciap=None, zfinger=None), TRAF_0)
    Initial(cIAP(btraf=None), cIAP_0)
    Initial(NSC(bnfkb=None), NSC_0)
    Initial(NFkB(bnsc=None, state = 'I'), NFkB_0)
    Initial(CYLD(btraf=None, state = 'U'),CYLD_0)
    Initial(FADD(bDD=None, bDED1=None, bDED2=None), FADD_0)

def TNFa_to_ComplexI():
    """Reaction network that produces Complex I and activates NFkB [1]. Recruitment of LUBAC and NFkB pathway components to complex I is approximated by a one-step conversion to "NFkB Signaling Complex" (reaction 7). A20 ubiquitylates RIP1 and targets it for degradation. The action of A20 is apporixmated at a single reaction (reaction 8). 
        
        1. TNFa + TNFR1 <> TNFa:TNFR1
        2. TNFa:TNFR1 + TRADD <> TNFa:TNFR1:TRADD
        3. TNFa:TNFR1:TRADD + RIP1(unmod) <> TNFa:TNFR1:TRADD:RIP1(unmod)
        4. TNFa:TNFR1:TRADD:RIP1(unmod) + TRAF2 <> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2
        5. TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2 + cIAP <> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:cIAP
        
        6. TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:cIAP >>TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:cIAP
        7. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 >> NFkB Signaling Complex
        8. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 >> TNFa:TNFR1:TRADD + TRAF2 
        9. TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2 + CYLD <> TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:CYLD
        10.TNFa:TNFR1:TRADD:RIP1(poly-ubiquitylated):TRAF2:CYLD >> TNFa:TNFR1:TRADD:RIP1(unmod):TRAF2:CYLD
        
        1. Olivier Micheau, Jurg Tschopp, Induction of TNF Receptor I-Mediated Apoptosis via Two Sequential Signaling Complexes, Cell, Volume 114, Issue 2, 25 July 2003, Pages 181-190
    """
    bind(TNFa(brec = None),'brec', TNFR1(blig = None), 'blig', [1e-6, 1e-3])
    bind(TNFR1(blig = ANY, bDD = None), 'bDD', TRADD(bDD1 = None, bDD2 = None), 'bDD1', [1e-6, 1e-3])
    #we need to get bind_complex working! 
    Rule('RIP1_to_complex1', TNFR1(blig = ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = None) + RIP1(bDD = None, btraf = None, state =  'unmod') <> TNFR1(blig = ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = 1)%RIP1(bDD = 1, btraf = None, state =  'unmod'), Parameter('k1', 1e-6),Parameter('k2',1e-3))
    Rule('TRAF_to_complex1', TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY, btraf = None, state =  'unmod') + TRAF(brip=None) <> TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY, btraf = 1, state =  'unmod')%TRAF(brip=1), Parameter('k3', 1e-6), Parameter('k4',1e-3))
    bind(TRAF(bciap =  None), 'bciap', cIAP(btraf = None), 'btraf', [1e-6, 1e-3])
    bind(TRAF(zfinger =  None), 'zfinger', CYLD(btraf = None), 'btraf', [1e-6, 1e-3])
    
    ComplexI = TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)
    
    Rule('RIP1_Ubiquitination', ComplexI %RIP1(bDD = ANY, btraf = ANY, state =  'unmod')%TRAF(brip = ANY, bciap = ANY) >> ComplexI %RIP1(bDD = ANY, btraf = None, state =  'ub')%TRAF(brip = ANY, bciap = ANY), Parameter('RIP1_ubiq_kc', 1e-1))
    Rule('RIP1_Deubiquitination', ComplexI%RIP1(bDD = ANY, btraf = ANY, state =  'ub')%TRAF(brip = ANY, zfinger = ANY) >> ComplexI%RIP1(bDD = ANY, btraf = None, state =  'unmod')%TRAF(brip = ANY, zfinger = ANY), Parameter('RIP1_deub_kc', 1e-1))
    Rule('Establish_NFkB_Signaling_Complex', ComplexI%RIP1(bDD = ANY, btraf = ANY, state =  'ub')%TRAF(brip = ANY, bciap = ANY) >> NSC(bnfkb=None), Parameter('NSC_esta_kc', 1e-8))
    Rule('RIP1_Degradation', ComplexI%RIP1(bDD = ANY, btraf = ANY, state =  'ub')%TRAF(brip = ANY) >> ComplexI + TRAF(brip = None), Parameter('RIP1_degr_kc', 1e-1))

def CompI_TRADD_RIP1_Dissociation():
    """Some believe that TRADD, RIP1 and others are released into the cytoplasm after RIP1 deubiquitination or after TNFR endocytosis. This mechanism presents a problem: The RIP1 released into the cytoplasm is in theory the same as RIP1 that existed in the cytoplasm prior to TNF stimulation. Meaning apoptosis whould spontaneously occur. Dissociation of TRADD:RIP1 from complex I is a hypothetical pathway that could distinguish pre- and post- TNF RIP1. It could also serve as a mechanism for FADD independent necrosome formation. This hypothesis is supported by a suggested mechanism, to explain signal transduction from TNFa to Riptosome formation[1]. 
        
        TNFa:TNFR1:TRADD:RIP1(unmod) <> TNFa:TNFR1 + TRADD:RIP1(unmod)
        
        1. Dickens, LS., IR Powley, MA Hughes, M MacFarlane, Exp. Cell. Res. 318 (2012) 1269-1277"""
    
    Rule('RIP1_TRADD_Complex1', TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = None) + TRADD(bDD1 = None, bDD2 = ANY)%RIP1(bDD = ANY, state =  'unmod') <> TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = 1)%TRADD(bDD1 = 1, bDD2 = ANY)%RIP1(bDD = ANY, state =  'unmod'), Parameter('k5', 1e-6),Parameter('k6',1e-3))

def CompII_Hypothesis_1():
    """Michaeu and Tschopp, showed that FADD is absent from Complex I[1]. The presence of TRADD in Complex II is debated [2]. As explained in CompI_TRADD_RIP1_Dissociation(), simply releasing RIP1 from Complex I so that it can bind FADD creates a condition where cell death spontaneously occurs. Here, I hypothesize that FADD transiently associates with Complex I in order to retrieve RIP1 from Complex I. This association, being transient, would have been difficult to detect.
        
        1. Olivier Micheau, Jurg Tschopp, Induction of TNF Receptor I-Mediated Apoptosis via Two Sequential Signaling Complexes, Cell, Volume 114, Issue 2, 25 July 2003, Pages 181-190
        
        2. Dickens, LS., IR Powley, MA Hughes, M MacFarlane, "The 'complexities' of life and death: Death receptor signalling platforms" Exp. Cell. Res. 318 (2012) 1269-1277"""

    Rule('TNFR1_to_FADD_1', TNFR1(bDD = 2)%TRADD(bDD1 = 2, bDD2 = 3)%RIP1(bDD = 3, btraf = 4)%TRAF(brip = 4) + FADD(bDD = None)>> TNFR1(bDD = 2)%TRADD(bDD1 = 2, bDD2 = None) + RIP1(bDD = 1, btraf = None)%FADD(bDD = 1) + TRAF(brip = None), Parameter('TNFR1_FADD_kc_1', 1e-1))
    
    Rule('TNFR1_to_FADD_2', TNFR1(bDD = 2)%TRADD(bDD1 = 2, bDD2 = 3)%RIP1(bDD = 3, btraf = None)+ FADD(bDD = None)>> TNFR1(bDD = 2)%TRADD(bDD1 = 2, bDD2 = None) + RIP1(bDD = 1, btraf = None)%FADD(bDD = 1), Parameter('TNFR1_FADD_kc_2', 1e-1))

def CompII_Hypothesis_2():
    """Dickens et al, speculated that TRADD:RIP1 is released into the cytoplasm following RIP1 deubiquitination by CYLD[1]. Furter, Dinckens question presence of TRADD in Complex II[1], but it is established that Complex II is contains FADD. To model these findings, we suggest that FADD replaces the TRADD in the cytoplasmic TRADD:RIP1 complex. """
        
    Rule('TRADD_to_FADD_1', TRADD(bDD1 = None, bDD2 = 2)%RIP1(bDD = 2)%TRAF(brip = 3) + FADD(bDD = None)>> TRADD(bDD1 = None, bDD2 = None) + RIP1(bDD = 1)%FADD(bDD = 1) + TRAF(brip = None), Parameter('TRADD_FADD_kc_1', 1e-1))
    
    Rule('TRADD_to_FADD_2', TRADD(bDD1 = None, bDD2 = 2)%RIP1(bDD = 2, btraf = None)+ FADD(bDD = None)>> TRADD(bDD1 = None, bDD2 = None) + RIP1(bDD = 1, btraf = None)%FADD(bDD = 1), Parameter('TRADD_FADD_kc_2', 1e-1))

def CompII_Hypothesis_3():
    """Dickens et al, speculated that TRADD:RIP1 is released into the cytoplasm following RIP1 deubiquitination by CYLD[1]. FADD is recruited to RIP1 to form a cytoplasmic complex II. The presence of TRADD in complex II is debated. Here, we hypothesize that it is present. """
    
    Rule('TRADD_to_FADD_3', FADD(bDD = None) + TRADD(bDD1 = None, bDD2 = ANY)%RIP1(bDD = ANY, btraf = None) <> FADD(bDD = 1)%TRADD(bDD1 = 1, bDD2 = ANY)%RIP1(bDD = ANY, btraf = None), Parameter('TRADD_FADD_kf', 1e-1), Parameter('TRADD_FADD_kr', 1e-1))
    
def ComplexII_to_Bid_Monomers():
    """ Bid is declared in Albeck Modules"""
    Monomer('proC8', ['bDED'])
    Monomer('C8', ['bf', 'state'], {'state':['A', 'I']})        #active caspase 8
    Monomer('flip_L',['bDED'])
    Monomer('flip_S',['bDED'])
    Monomer('RIP3', ['bRHIM', 'state'], {'state':['unmod', 'po4', 'trunc', 'N']})
    Monomer('MLKL', ['bRHIM', 'state'], {'state':['unmod', 'active', 'inactive']})

def ComplexII_to_Bid_Initials():
    Parameter('flip_L_0',  39022) # 39022 molecules per cell (J Immunol 2005)
    Parameter('flip_S_0',  39022) # molecules per cell assummed both isoforms w/conc.
    Parameter('proC8_0' ,  16057) # procaspase 8 molecules per cell 16057 (J Immunol 2005)
    Parameter('C8_0'    ,      0) # active caspase 8 dimers per cell.
    Parameter('RIP3_0'  ,  2.0e4) # molecules per cell
    Parameter('MLKL_0'  , 1.0e6) # molecules per cell
    
    alias_model_components()
    Initial(RIP3(bRHIM = None, state = 'unmod'), RIP3_0)# RIP3
    Initial(flip_L(bDED = None), flip_L_0)
    Initial(flip_S(bDED = None), flip_S_0)
    Initial(proC8(bDED = None), proC8_0)
    Initial(C8(bf = None, state = 'I'), C8_0)
    Initial(MLKL(bRHIM = None, state = 'unmod'), MLKL_0)

def ComplexIIa_Assembly():
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
    #   X = RIP1, TRADD

    # ----------procaspase8 and cFlip recruitment-----
    bind(FADD(bDD = ANY, bDED1 = None),'bDED1', proC8(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    Rule('bind_FADD_proC8_2', FADD(bDD = ANY, bDED2 = None) + proC8(bDED = None) <> FADD(bDD = ANY, bDED2 = 1)%proC8(bDED = 1), Parameter('bind_FADD_proC8_2_kf', 7.27e-06), Parameter('bind_FADD_proC8_2_kr', 0.018)) # (J Immunol 2005) FIX: Bind function fails when one molecule can bind another at more than one binding site. The failure occurs because _rule_name_generic does not distinguish between binding sites.
    
    bind(FADD(bDD = ANY, bDED1 = None),'bDED1', flip_L(bDED = None), 'bDED', [7.27e-05, 0.018]) # (J Immunol 2005)
    Rule('bind_FADD_flip_L_2', FADD(bDD = ANY, bDED2 = None) + flip_L(bDED = None) <> FADD(bDD = ANY, bDED2 = 1)%flip_L(bDED = 1), Parameter('bind_FADD_flip_L_2_kf', 7.27e-06), Parameter('bind_FADD_flip_L_2_kr', 0.018)) # (J Immunol 2005)[...source] says cFlip-L has greater affinity for proC8 than proC8 itself.
    
    bind(FADD(bDD = ANY, bDED1 = None),'bDED1', flip_S(bDED = None), 'bDED', [7.27e-06, 0.018]) # (J Immunol 2005)
    Rule('bind_FADD_flip_S_2', FADD(bDD = ANY, bDED2 = None) + flip_S(bDED = None) <> FADD(bDD = ANY, bDED2 = 1)%flip_S(bDED = 1), Parameter('bind_FADD_flip_S_2_kf', 7.27e-06), Parameter('bind_FADD_flip_S_2_kr', 0.018)) # (J Immunol 2005)

def ComplexIIb_to_MLKL():
    """Necrosome formation and MLKL activation"""
    bind(RIP1(bDD = ANY, bRHIM = None, state = 'unmod'), 'bRHIM', RIP3(bRHIM = None, state = 'unmod'), 'bRHIM', [1e-6, 1e-3])
    
    Rule('Rip3_PO4lation', RIP1(bRHIM=ANY, state = 'unmod')%RIP3(bRHIM=ANY, state='unmod') >> RIP1(bRHIM=ANY, state = 'unmod')%RIP3(bRHIM=ANY, state = 'po4'), Parameter('k19', 1e-2))
    Rule('Rip1_PO4lation', RIP1(bRHIM=ANY, state = 'unmod')%RIP3(bRHIM=ANY, state='po4') >> RIP1(bRHIM=ANY, state = 'po4')%RIP3(bRHIM=ANY, state = 'po4'), Parameter('k20', 1e-3))
    
    catalyze_state(RIP1(state='po4'), 'bMLKL', MLKL(), 'bRHIM', 'state', 'unmod', 'active', [1e-6,1e-3, 1e-1])
    catalyze_state(MLKL(state='active'), 'bRHIM', MLKL(), 'bRHIM', 'state', 'unmod', 'active', [1e-7, 0.4, 0.01])

def RIP1_truncation_ComplexII():
    #Flip_L:proC8 mediated RIP1 truncation
    Rule('RIP1_trunc_1a', RIP1(bDD = None, state = 'unmod') + FADD(bDD = None, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED=ANY)%proC8(bDED=ANY) <> RIP1(bDD = 1, state = 'unmod')%FADD(bDD = 1, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED=ANY)%proC8(bDED=ANY), Parameter('RIP1_trunc_kcf1', 1e-6), Parameter('RIP1_trunc_kcr1', 1e-3))

    Rule('RIP1_trunc_1b', RIP1(bDD = 1, state = 'unmod')%FADD(bDD = 1, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED=ANY)%proC8(bDED=ANY) >> RIP1(bDD = None, state = 'trunc') + FADD(bDD = None, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED=ANY)%proC8(bDED=ANY), Parameter('RIP1_trunc_kc1', 1e-1))

    Rule('RIP1_trunc_2a', RIP1(bDD = None, state = 'unmod') + TRADD(bDD1 = 2, bDD2 = None)%FADD(bDD = 2, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED=ANY)%proC8(bDED=ANY) <> RIP1(bDD = 1, state = 'unmod')%TRADD(bDD1 = 2, bDD2 = 1)%FADD(bDD = 2, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED=ANY)%proC8(bDED=ANY), Parameter('RIP1_trunc_kcf2', 1e-6), Parameter('RIP1_trunc_kcr2', 1e-3))

    Rule('RIP1_trunc_2b', RIP1(bDD = ANY, state = 'unmod')%TRADD(bDD1 = ANY, bDD2 = ANY)%FADD(bDD = ANY, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED=ANY)%proC8(bDED=ANY) >> RIP1(bDD = None, state = 'trunc') + TRADD(bDD1 = ANY, bDD2 = None)%FADD(bDD = ANY, bDED1 = ANY, bDED2 = ANY)%flip_L(bDED=ANY)%proC8(bDED=ANY), Parameter('RIP1_trunc_kc2', 1e-1))

    #proC8:proC8 mediated RIP1 truncation
    Rule('RIP1_trunc_3a', RIP1(bDD = None, state = 'unmod') + FADD(bDD = None, bDED1 = 2, bDED2 = 3)%proC8(bDED=2)%proC8(bDED=3) <> RIP1(bDD = 1, state = 'unmod')%FADD(bDD = 1, bDED1 = 2, bDED2 = 3)%proC8(bDED=2)%proC8(bDED=3), Parameter('RIP1_trunc_kcf3', 1e-6), Parameter('RIP1_trunc_kcr3', 1e-3))

    Rule('RIP1_trunc_3b', RIP1(bDD = ANY, state = 'unmod')%FADD(bDD = ANY, bDED1 = ANY, bDED2 = ANY)%proC8(bDED=ANY)%proC8(bDED=ANY) >> RIP1(bDD = None, state = 'trunc') + FADD(bDD = None, bDED1 = None, bDED2 = None) + C8(bf = None, state = 'A'), Parameter('RIP1_trunc_kc3', 1e-1))

    Rule('RIP1_trunc_4a', RIP1(bDD = None, state = 'unmod') + TRADD(bDD1 = 2, bDD2 = None)%FADD(bDD = 2, bDED1 = 3, bDED2 = 4)%proC8(bDED=3)%proC8(bDED=4) <> RIP1(bDD = 1, state = 'unmod')%TRADD(bDD1 = 2, bDD2 = 1)%FADD(bDD = 2, bDED1 = 3, bDED2 = 4)%proC8(bDED=3)%proC8(bDED=4), Parameter('RIP1_trunc_kcf4', 1e-6), Parameter('RIP1_trunc_kcr4', 1e-3))

    Rule('RIP1_trunc_4b', RIP1(bDD = ANY, state = 'unmod')%TRADD(bDD1 = ANY, bDD2 = ANY)%FADD(bDD = ANY, bDED1 = ANY, bDED2 = ANY)%proC8(bDED=ANY)%proC8(bDED=ANY) >> RIP1(bDD = None, state = 'trunc') + TRADD(bDD1 = ANY, bDD2 = None)%FADD(bDD = ANY, bDED1 = None, bDED2 = None)+ C8(bf = None, state = 'A'), Parameter('RIP1_trunc_kc4', 1e-1))

    # Caspase 8 activation
    Rule('C8_activation1', FADD(bDD = None, bDED2 = ANY, bDED1 = ANY)%proC8(bDED = ANY)%proC8(bDED = ANY) >> FADD(bDD = None, bDED1 = None, bDED2 = None) + C8(bf = None, state = 'A'), Parameter('kc_c8_1', 1))

    Rule('C8_activation2', TRADD(bDD1 = None, bDD2 = 1)%FADD(bDD = 1, bDED2 = ANY, bDED1 = ANY)%proC8(bDED = ANY)%proC8(bDED = ANY) >> TRADD(bDD1 = None, bDD2 = 1)%FADD(bDD = 1, bDED1 = None, bDED2 = None) + C8(bf = None, state = 'A'), Parameter('kc_c8_2', 1))

def C8_catalyzed_truncations():
    
    #   RIP1 + C8 <-> RIP1:C8 >> RIP1[trunc] + C8
    #   RIP3 + C8 <-> RIP3:C8 >> RIP3[trunc] + C8
    #   Bid + C8 <-> Bid:C8 >> Bid[trunc] + C8
    #   CYLD + C8 <-> CYLD:C8 >> CYLD
    
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP1(bDD=None), 'bRHIM', 'state', 'unmod', 'trunc', [1e-6, 1e-3, 1e-1])
    
    #RIP3 Truncation
    catalyze_state(C8(bf = None, state = 'A'), 'bf', RIP3(), 'bRHIM', 'state', 'unmod', 'trunc', [1e-6, 1e-3, 1e-1])
    
    #Bid Truncation
    catalyze_state(C8(bf = None, state = 'A'), 'bf', Bid(), 'bf', 'state', 'U', 'T', [1.e-6, 0.001, 0.1])

    #CYLD Truncation
    catalyze_state(C8(bf = None, state = 'A'), 'bf', CYLD(), 'btraf', 'state', 'U', 'T', [1e-6, 1e-3, 0.1])

def NFkB_Activation_and_Signaling_monomers():
    Monomer('Flip_L_degradase', ['bf'])
    Monomer('Flip_S_degradase', ['bf'])
    Monomer('TRAF_degradase', ['bf'])


def NFkB_Activation_and_Signaling_Initials():
    """Several papers mention NFkB mediated expression of c-Flip as an important mediator of apoptosis and necrosis."""
    Parameter('Flip_L_degradase_0', 0)
    Parameter('Flip_S_degradase_0', 0)
    Parameter('TRAF_degradase_0', 0)
    alias_model_components()
    
    Initial(Flip_L_degradase(bf=None), Flip_L_degradase_0)
    Initial(Flip_S_degradase(bf=None), Flip_S_degradase_0)
    Initial(TRAF_degradase(bf=None), TRAF_degradase_0)

def NFkB_Activation_and_Signaling():
    """This model focuses on the decision between apoptotic and necroptotic cell death pathways. The NFkB pathway, roughly approximated here, does not induce cell death. But it may influence the expression of key regulators of the cell death pathways."""
    catalyze_state(NSC(bnfkb=None), 'bnfkb', NFkB(state = 'I'), 'bnsc', 'state', 'I', 'A', [1e-6, 1e-3, 0.1])
    
    Rule('NFkB_cFlipL', NFkB(state = 'A') >> NFkB(state = 'A') + flip_L(bDED=None), Parameter('NFkB_FlipL', 1e-2))
    Rule('NFkB_cFlipS', NFkB(state = 'A') >> NFkB(state = 'A') + flip_S(bDED=None), Parameter('NFkB_FlipS', 1e-2))
    Rule('NFkB_TRAF', NFkB(state = 'A') >> NFkB(state = 'A') + TRAF(brip=None, bciap=None, zfinger=None), Parameter('NFkB_TRAF_kc', 1e-2))
    
    Rule('NFkB_FlipL_deg', NFkB(state = 'A') >> NFkB(state = 'A') + Flip_L_degradase(bf=None), Parameter('Degradase_flipL', 1e-6))
    Rule('NFkB_FlipS_deg', NFkB(state = 'A') >> NFkB(state = 'A') + Flip_S_degradase(bf=None), Parameter('Degradase_flipS', 1e-6))
    Rule('NFkB_TRAF_deg', NFkB(state = 'A') >> NFkB(state = 'A') + TRAF_degradase(bf=None), Parameter('Degradase_TRAF', 1e-6))
    
    Rule('Deg_cFlipL', Flip_L_degradase(bf=None) + flip_L(bDED=None) >> Flip_L_degradase(bf=None), Parameter('deg_FlipL', 5e-6))
    Rule('Deg_cFlipS', Flip_S_degradase(bf=None) + flip_S(bDED=None) >> Flip_S_degradase(bf=None), Parameter('deg_FlipS', 5e-6))
    Rule('Deg_TRAF', TRAF_degradase(bf=None) + TRAF(brip=None, bciap=None, zfinger=None) >> TRAF_degradase(bf=None), Parameter('deg_TRAF', 5e-6))

def Bid_Hypothesis_monomers():
    Monomer('BidK', ['bf']) #unknown Bid-kinase

def Bid_Hypothesis_initials():
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
    
    Parameter('BidK_0'  , 5.0e3) # molecules per cell
    
    alias_model_components()
    Initial(BidK(bf = None), BidK_0)

def Bid_Hypothesis():
    """The Zinkel lab found evidence of Bid mediated inhibition of necroptosis. Dr. Zinkel proposed that Bid has a third state (e.g. phosporylated) that can inhibit the necrosome by sequestering RIP3. We presume that Bid-po4 cannot be truncated. There are 3 phosphorylation sites on Bid and evidence that BId-po4 is resistant to truncation exists"""
    
    #-------------Bid Interactions--------------------------
    # Bid Phosphorylation and Truncation
    catalyze_state(BidK(), 'bf', Bid(), 'bf', 'state', 'U', 'po4', [1e-6, 1e-3, 1e-1])
    
    # Bid-PO4 sequestering RIP1
    bind(RIP1(bDD = None, bRHIM = None, state = 'unmod'), 'bRHIM', Bid(bf = None, state = 'po4'), 'bf', [1e-6, 1e-3])

def Bid_RIP_recruits_proC8():
    #Bid mediated inhibition of necrosis (hypothesis 1)
    Rule('Bid_recruits_1', RIP1(bDD = None, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') + proC8(bDED = None) <> RIP1(bDD = 1, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') % proC8(bDED = 1), Parameter('kbid1', 1e-6), Parameter('kbid2', 1e-3))

def Bid_RIP_proC8_truncation():
    Rule('Bid_trunc_RIP_1', RIP1(bDD = ANY, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') % proC8(bDED = ANY) >> Bid(bf = None, state = 'po4') + RIP1(bDD = None, bRHIM = None, state = 'trunc') + proC8(bDED = None), Parameter('kbid5', 1e-1))

def Bid_RIP_proC8_to_NFkB():
    #Bid mediated inhibition of necrosis (hypothesis 1)
    Rule('Bid_RIP_NFkB', RIP1(bDD = ANY, bRHIM =  ANY, state = 'unmod') % Bid(bf = ANY, state = 'po4') % proC8(bDED = ANY) >> NSC(bnfkb=None), Parameter('kbid7', 1e-1))

def C3_inhibits_MLKL():
    """it has not been established whether MLKL is de-activated in apoptosis. But inorder to make a cellular decision you have to have one choice, once selected, inhibit all of the alternative pathways. This is a hypothetical apoptosis mediated inhibition of necrosis effector, MLKL."""
    catalyze_state(C3(state='A'), 'bf', MLKL(), 'bRHIM', 'state','unmod','inactive', [1e-6, 1e-3, 1])

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
    Observable('Obs_TNFa', TNFa(brec =  None))
    Observable('ComplexI', TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(state = 'unmod'))
    Observable('ComplexI_ub', TNFa(brec = ANY)%TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(state = 'ub'))
    Observable('ComplexI_TRAF', TNFR1(blig =  ANY, bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY, btraf = 1, state =  'unmod')%TRAF(brip=1))
    Observable('TRADD_RIP1', TRADD(bDD1 = None, bDD2 = ANY)%RIP1(bDD = ANY, state =  'unmod'))
    Observable('TRADD_RIP1_2', TRADD(bDD1 = None, bDD2 = ANY)%RIP1(bDD = ANY, btraf = None))
    Observable('ComplexII', FADD(bDD = ANY)%TRADD(bDD1 = ANY, bDD2 = ANY)%RIP1(bDD = ANY, btraf = None))
    Observable('ComplexIIa', FADD(bDD = ANY, bDED1 = ANY, bDED2=ANY))
    Observable('Obs_NFkB', NFkB(state = 'A'))
    Observable('Obs_FADD_Sole', FADD(bDD = None))
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
