######################################################################
#
# Lilith routine example
# To put in Lilith-2.X/examples/python/ folder 
# To execute from /Lilith-2.X root folder
#
# Compute the SM point compatibility in the (CGa, Cg) model
#
# Use the library iminuit to perform minimization
# Use the library scipy to compute p-values
#
######################################################################

import sys, os
import iminuit
from iminuit import Minuit
from math import floor
from scipy import stats

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

iminuit_version = iminuit.__version__
print("iminuit version:", iminuit_version)
iminuit_version_f = float(iminuit_version[:3])

######################################################################
# Parameters
######################################################################

# Choose which experimental data to fit
#myexpinput = "data/latest.list"
myexpinput = "data/latestRun2.list"
# Lilith precision mode
myprecision = "BEST-QCD"
# Higgs mass to test
myhmass = 125.09

verbose=False
timer=False

print("\nTask: Compute the SM point compatibility in the (CGa, Cg) model")
print("Experimental input:", myexpinput)

######################################################################
# * usrXMLinput: generate XML user input
# * getL :       -2LogL for a given (CGa, Cg)
######################################################################

def usrXMLinput(mass=125., CGa=1, Cg=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings
       CGa, Cg"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="ff">1.</C>
  <C to="VV">1.</C>
  <C to="gammagamma">%(CGa)s</C>
  <C to="gg">%(Cg)s</C>

  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass':mass, 'CGa':CGa, 'Cg':Cg,'precision':precision}
        
    return myInputTemplate%myInput



def getL_CGaCg(CGa, Cg):
  """ 
      Returns -2Log(L) for a model in which BSM loop-contributions 
      to the gamma-gamma and gluon-gluon processes are present. 
      The tree-level couplings CU, CD, CV are assumed to be 1.
  """
  myXML_user_input = usrXMLinput(myhmass, CGa=CGa, Cg=Cg, precision=myprecision)
  lilithcalc.computelikelihood(userinput=myXML_user_input)
  return lilithcalc.l


######################################################################
# Calculations
######################################################################

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose,timer)
# Read experimental data
lilithcalc.readexpinput(myexpinput)

print("\n***** evaluating SM -2LogL *****")
# Evaluate Likelihood at the SM point
SM_minus2logL = getL_CGaCg(1., 1.)
print("-2LogL(SM) =", SM_minus2logL)

print("\n***** performing (CGa,Cg) model fit *****")
# Initialize the fit; parameter starting values and limits
if iminuit_version_f < 2.0:
  m = Minuit(getL_CGaCg, CGa=0.9, limit_CGa=(0,3), Cg=0.9, limit_Cg=(0,3), print_level=0, errordef=1, error_CGa=1, error_Cg=1)
else:
  m = Minuit(getL_CGaCg, CGa=0.9, Cg=0.9)
  m.limits = [(0, 3), (0, 3)]
  m.errordef = 1 # 1 for -2LogL (or least square), 0.5 for -LogL 
  m.errors = [1, 1]
  m.print_level = 0

# Fit the model
m.migrad()
# Display parameter values at the best-fit point
print("Best-fit point: ")
print("Cgamma =", m.values["CGa"])  
print("Cgluon =", m.values["Cg"]) 

bestfit_CGa_Cg_minus2logL = m.fval

print("-2LogL =", bestfit_CGa_Cg_minus2logL)

print("\n***** model comparison *****")

delta_minus2logL = SM_minus2logL - bestfit_CGa_Cg_minus2logL
pval = 1 - stats.chi2.cdf(delta_minus2logL, 2)

print("Delta(-2logL) =", delta_minus2logL)
print("p-value =", pval, "\n")


