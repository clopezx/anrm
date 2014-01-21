from pysb import *
from pysb.util import alias_model_components
from earm import shared

from anrm import irvin_modules_v2 as irvin

Model('model')

irvin.Momomers_FasL_to_DISC()
irvin.Initials_Fas_to_DISC()
irvin.CD95_to_SecondaryComplex()
irvin.FADD_to_C8()

"""irvin.TNFa_to_ComplexI_Monomers()
irvin.TNFa_to_ComplexI_Initials()
irvin.TNFa_to_ComplexI()
irvin.CompI_TRADD_RIP1_Dissociation()
#irvin.CompII_Hypothesis_1() #FADD transiently localizes to TNFR1 to retrieve RIP1"""