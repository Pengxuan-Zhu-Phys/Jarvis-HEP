###############################################################
#
# Lilith routine example
# To put in Lilith-2.X/examples/python/ folder 
# To execute from /Lilith-2.X root folder
#
# Constraints on reduced couplings in the (CV, CF) model. 
# CV: common coupling to EW gauge bosons,
# CF: common coupling for fermions.
# The loop-induced effective couplings CGa and Cg are
# determined from CV and CF.
#
# 1-dimensional likelihood profiles are obtained from a
# profile-likelihood analysis
#
# Use the libraries matplotlib, iminuit
#
###############################################################

import sys, os
import matplotlib.pyplot as plt
import iminuit
from iminuit import Minuit
import matplotlib

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

# Exprimental results
#myexpinput = "data/latest.list"
myexpinput = "data/latestRun2.list"

# Lilith precision mode
myprecision = "BEST-QCD"

# Output
if (not os.path.exists("results")):
    os.mkdir("results")
outputplot = "results/CVCF_1dprofiles.pdf"

# Higgs mass to test
myhmass = 125.09

verbose=False
timer=False

print("\nTask: profile-likelihood analysis in (CF,CV) model")
print("Experimental input:", myexpinput)

######################################################################
# * usrXMLinput: generate XML user input
# * getL:        -2LogL for a given (CV, CF) point
######################################################################

def usrXMLinput(mh=125., CV=1., CF=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings
       the gamma-gamma and/or gluon-gluon couplings are 
       evaluated as functions of the fermionic and bosonic
       reduced couplings """
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mh)s</mass>
  
  <C to="tt">%(CF)s</C>
  <C to="bb">%(CF)s</C>
  <C to="tautau">%(CF)s</C>
  <C to="cc">%(CF)s</C>
  <C to="WW">%(CV)s</C>
  <C to="ZZ">%(CV)s</C>
  
  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>
  
  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mh':mh, 'CV':CV, 'CF':CF, 'precision':precision}
        
    return myInputTemplate%myInput



def getL(CV, CF):
    myXML_user_input = usrXMLinput(mh=myhmass, CV=CV, CF=CF, precision=myprecision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    return lilithcalc.l


######################################################################
# Calculations
######################################################################

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose, timer)
# Read experimental data
lilithcalc.readexpinput(myexpinput)

print("\n***** performing model fit: migrad, hesse, minos *****")

# Initialize the fit; parameter starting values and limits
if iminuit_version_f < 2.0:
  m = Minuit(getL, CV=1, limit_CV=(0,3), CF=1, limit_CF=(0,3), print_level=0, errordef=1, error_CV=0.2, error_CF=0.2)
else:
  m = Minuit(getL, CV=1, CF=1)
  m.limits = [(0, 3), (0, 3)]
  m.errordef = 1 # 1 for -2LogL (or least square), 0.5 for -LogL 
  m.errors = [0.2, 0.2]
  m.print_level = 0

# Minimization and error estimation
m.migrad()
m.hesse()   # run covariance estimator
m.minos()


print("\nbest-fit point:", m.values) 
print("\nHesse errors:", m.errors)
print("\nMinos errors:")
for key, value in list(m.merrors.items()):
    print(key, value)

if iminuit_version_f < 2.0:
  print("\nCorrelation matrix:\n", m.matrix(correlation=True))
  print("\nCovariance matrix:\n", m.matrix())
else:
  print("\nCorrelation matrix:\n", m.covariance.correlation())
  print("\nCovariance matrix:\n", m.covariance)

# Display parameter values at the best-fit point
#print "\nBest-fit point:" 
#print "CV =", m.values["CV"], "\nCF =", m.values["CF"], "\n-2LogL =", m.fval


print("\n***** getting the 1-dimensional likelihood profile of CV *****")
# Profiling over CF
if iminuit_version_f < 2.0:
  xV,yV,rV = m.mnprofile('CV', bins=300, bound=(0., 3), subtract_min=True)
else:
  xV,yV,rV = m.mnprofile('CV', size=300, bound=(0., 3), subtract_min=True)

print("***** getting the 1-dimensional likelihood profile of CF *****")
# Profiling over CV
if iminuit_version_f < 2.0:
  xF,yF,rF = m.mnprofile('CF', bins=300, bound=(0., 3), subtract_min=True)
else:
  xF,yF,rF = m.mnprofile('CF', size=300, bound=(0., 3), subtract_min=True)


######################################################################
# Plot
######################################################################

print("***** plotting *****")

matplotlib.rcParams['xtick.major.pad'] = 18
matplotlib.rcParams['ytick.major.pad'] = 18

fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(121)

plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xV,yV,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([0.8,1.3])
plt.xlabel(r'$C_V$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

plt.title("Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=20, ha="left")

ax = fig.add_subplot(122)

plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xF,yF,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([0.7,1.4])
plt.xlabel(r'$C_F$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

plt.title("Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=20, ha="left")

#fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("result see " + lilith_dir + "/" + outputplot)
print("***** done *****\n")



