###############################################################
#
# Lilith routine example
# To execute from /Lilith-2.X root folder
#
# Constraints on stau masses and mixing angle
# from the additionnal loop contribution to
# the Higgs -> gamma gamma decay rate
# (exemple from the  Lilith manual arXiv:1502.04138)
#
# Use the libraries matplotlib and numpy
#
###############################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from math import floor, log, sqrt, sin, cos, atan
from cmath import sqrt as csqrt
from cmath import asin as casin

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith


######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = "data/latestRun2.list"
# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# Stau 1 physical mass [GeV]
m1_fixed = 85
# Tan beta
tb = 10

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/grid_staus_theta_m2.out"
outputplot = "results/m1stau"+str(m1_fixed)+"GeV_tanbeta"\
              +str(tb)+"_theta_m2stau.pdf"

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100


########################################################################
# SM and staus LO contribution to the reduced H-gamma-gamma coupling CGa
# Tree-level reduced couplings are assumed to be equal to 1
########################################################################

pi = np.pi
mW = 80.398
mZ = 91.1876

mt = 173.1
mb = 4.75
mc = 1.4
mtau = 1.777
sW2 = 0.23116

def fhiggs(t):
    if t<=1.:
        return casin(sqrt(t))**2.
    else:
        return -(log((sqrt(t)+sqrt(t-1.))/(sqrt(t)-sqrt(t-1.)))-pi*1j)**2./4.

def A0(tau):
    return -1./tau *(1.-1./tau * fhiggs(tau))

def A12(tau):
    return 2./tau *(1.+(1.-1./tau) * fhiggs(tau))

def A1(tau):
    return -(3.*tau+2.*tau**2. +3.*(2.*tau-1.) * fhiggs(tau))/tau**2


def get_CGa(m1, m2, theta):
  """ 
      Returns CGa computed from the SM particles contribution plus 
      staus contribution assuming decoupling of other states.
  """

  A12t = A12((hmass/(2.*mt))**2.)
  A12c = A12((hmass/(2.*mc))**2.)
  A12b = A12((hmass/(2.*mb))**2.)
  A12tau = A12((hmass/(2.*mtau))**2.)
  A1W = A1((hmass/(2.*mW))**2.)

  DL = cos(2*atan(tb))*(-1/2.+sW2)*mZ**2
  DR = cos(2*atan(tb))*(-sW2)*mZ**2
  d_mLL2 = (mtau**2. + DL)
  d_mRR2 = (mtau**2. + DR)
  d_mLR2 = 1/2.*1/2.*(m1**2-m2**2)*sin(2*theta)
  
  g1 = 1/2.*(d_mLL2+d_mRR2) + 1/2.*(d_mLL2-d_mRR2)*cos(2*theta) + d_mLR2*sin(2*theta)
  g2 = 1/2.*(d_mLL2+d_mRR2) - 1/2.*(d_mLL2-d_mRR2)*cos(2*theta) - d_mLR2*sin(2*theta)

  SM_amplitude = 4./3.*(A12t + A12c) + 1./3.*A12b + A12tau + A1W
  stau_amplitude = g1/m1**2.*A0((hmass/(2.*m1))**2.) \
                 + g2/m2**2.*A0((hmass/(2.*m2))**2.)
  
  return sqrt( abs(SM_amplitude + stau_amplitude)**2./abs(SM_amplitude)**2. )



######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125., CGa=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings
       CGa, Cg"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="ff">1.</C>
  <C to="VV">1.</C>
  <C to="gammagamma">%(CGa)s</C>
  <C to="gg">1.</C>

  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass':mass, 'CGa':CGa,'precision':precision}
        
    return myInputTemplate%myInput

######################################################################
# Scan initialization
######################################################################

print("***** stau1 mass =", str(m1_fixed), "GeV, tan(beta) =", str(tb), "*****\n")
print("***** scan initialization *****")

m2_min = max(100, m1_fixed)
m2_max = 1000.
theta_min = 0.
theta_max = pi/2.

# Prepare output
fresults = open(output, 'w')

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose=False,timer=False)
# Read experimental data
lilithcalc.readexpinput(exp_input)


######################################################################
# Scan routine
######################################################################

m2logLmin=10000
max=-1

print("***** running scan *****")

for theta in np.linspace(theta_min, theta_max, grid_subdivisions):
    fresults.write('\n')
    for m2 in np.linspace(m2_min, m2_max, grid_subdivisions):
        CGa = get_CGa(m1_fixed, m2, theta)
        myXML_user_input = usrXMLinput(hmass, CGa=CGa, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            m2min = m2
            thetamin = theta
        fresults.write('%.5f    '%theta +'%.5f    '%m2 + '%.5f     '%m2logL + '\n')

fresults.close()

print("***** scan finalized *****")

######################################################################
# Plot routine
######################################################################


print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig, ax = plt.subplots()

plt.minorticks_on()
plt.tick_params(labelsize=14, length=10, width=1.3)
plt.tick_params(which='minor', length=7, width=0.6)

# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,2]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
try:
# griddata using Natural Neighbor (nn) interpolation
    Z = griddata((x, y), z2, (X,Y), method='nearest')
except:
# If you don't have natgrid installed, use instead linear interpolation
    Z = griddata((x, y), z2, (X,Y), method="linear")

# heat map of -2LogL
cax = ax.imshow(Z, vmin=0, vmax=15, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], \
              aspect=(theta_max-theta_min)/(m2_max-m2_min), cmap=plt.get_cmap("jet_r"))

# Plotting the 68%, 95% and 99.7% CL contours
CS = ax.contour(xi,yi,Z,[2.3,5.99,11.83],linewidths=[2.5,2.5,2.5],colors=["#B22222","#FF8C00","#FFD700"])
ax.clabel(CS, inline=1, fontsize=11, colors='white')


# Title, labels, color bar...
plt.title("      Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=12, ha="left")
plt.text(1.3*pi/8., 184, r"$m_{\widetilde{\tau}_1}="+str(m1_fixed)+"$"r"$\rm{\ GeV}$, $\tan\beta="+str(tb)+"$", fontsize=16, backgroundcolor="white")
plt.xlabel(r'$\theta_{\widetilde{\tau}}$',fontsize=22)
plt.ylabel(r'$m_{\widetilde{\tau}_2} \rm{[GeV]}$',fontsize=22)
#plt.xticks([0, pi/8., pi/4., 3*pi/8., pi/2.],
#          ['$0$', r'$\frac{\pi}{8}$', r'$\frac{\pi}{4}$', r'$\frac{3\pi}{8}$', r'$\frac{\pi}{2}$'], fontsize=20)
plt.xticks([0, pi/8., pi/4., 3*pi/8., pi/2.],
          ['$0$', r'$\pi/8$', r'$\pi/4$', r'$3\pi/8$', r'$\pi/2$'], fontsize=18)

cbar = fig.colorbar(cax)
cbar.set_label(r"$\Delta(-2\log L)$",fontsize=20)

#plt.tight_layout()
fig.set_tight_layout(True)

# Saving figure (.pdf)
plt.savefig(outputplot)

print("results are stored in", lilith_dir + "/results")
print("***** done *****")
