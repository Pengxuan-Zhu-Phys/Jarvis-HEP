##################################################################
#
# Lilith routine example
#
# To put in Lilith-2.X/examples/python/ folder 
# To execute from /Lilith-2.X root folder
#
# Constraints on cos(beta-alpha) and tan(beta) in the 2HDM of
# Types I and II (in the sin(beta-alpha)>0 convention).
#
# Constraints are obtained through a likelihood-ratio test.
#
# Assumptions:
# * The lighter CP-even state h is identified with the observed
#   one (this fixes the coupling structure of the model)
#
# * The charged Higgs states are decoupled and thus
#   do not appear in the Higgs-photon-photon loop.
#   Only Higgs-fermion-fermion coupling modifications affect
#   this loop.
#
# Uses the libraries matplotlib (plotting) and numpy (functions)
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith


######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = "data/latest.list"
# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# 2HDM type = 1, 2
type = 1

# Number of steps for the square grid (cba,tb), need ~300 for fine grid but quite long to run
grid_subdivisions = 150

######################################################################

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")

if type == 1:
  output = "results/cba_tb_I_h_2d.out"
  outputplot = "results/cba_tb_I_h_2d.pdf"

if type == 2:
  output = "results/cba_tb_II_h_2d.out"
  outputplot = "results/cba_tb_II_h_2d.pdf"

# Scan ranges
if type == 1:
  cba_min = -0.5
  cba_max = 0.5
  tb_min = 0.1
  tb_max = 10.

if type == 2:
  cba_min = -0.1
  cba_max = 0.6
  tb_min = 0.1
  tb_max = 10.


######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, cba=0., tb=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings CU, CD, CV"""
    
    sba = np.sqrt(1-cba**2)
    
    if type == 1:
      CV = sba
      CU = sba + cba/tb
      CD = sba + cba/tb
    
    elif type == 2:
      CV = sba
      CU = sba + cba/tb
      CD = sba - cba*tb

    else:
      print("Error: 2HDM type parameter should be 1 or 2")
      sys.exit()

    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="uu">%(CU)s</C>
  <C to="dd">%(CD)s</C>
  <C to="VV">%(CV)s</C>
  
  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""

    myInput = {'mass':mass, 'CV':CV, 'CU':CU, 'CD':CD, 'precision':precision}
        
    return myInputTemplate%myInput

######################################################################
# Scan initialization
######################################################################

print("***** 2HDM scan initialization *****")

## Prepare output
fresults = open(output, 'w')
#
## Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose=False,timer=False)
## Read experimental data
lilithcalc.readexpinput(exp_input)
#
#
#######################################################################
## Scan routine
#######################################################################
#
m2logLmin=10000
max=-1

print("***** running 2HDM scan *****")

for cba in np.linspace(cba_min, cba_max, grid_subdivisions):
    fresults.write('\n')
    for tb in np.linspace(tb_min, tb_max, grid_subdivisions):
        myXML_user_input = usrXMLinput(hmass, tb=tb, cba=cba, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l                 #This is -2*Ln(L) at the (cba,tb) point
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            cbamin = cba
            tbmin = tb
        fresults.write('%.5f    '%cba +'%.5f    '%tb + '%.5f     '%m2logL + '\n')
fresults.close()

print("***** scan finalized *****")

######################################################################
# Plot routine
######################################################################


print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig = plt.figure(figsize=(7.5,7))
ax = fig.add_subplot(111)

plt.minorticks_on()
plt.tick_params(labelsize=20, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,2]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[el - z.min() for el in z]

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)

try:
# griddata using Natural Neighbor (nn) interpolation
    Z = griddata((x, y), z2, (X, Y), method="nearest")

except:
# If you don't have natgrid installed, use instead linear interpolation
    Z = griddata((x, y), z2, (X, Y), method="linear")

levels = np.arange(-2.0, 1.601, 0.4)
cmap = matplotlib.cm.YlOrRd

# Plotting the 68%, 95% and 99.7% CL contours (filled regions)
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],cmap=matplotlib.cm.get_cmap(cmap, len(levels) ))

#ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
#              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

# Title, labels, color bar etc.
plt.title("          Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=15, ha="left")
plt.xlabel(r'$\cos(\beta-\alpha)$',fontsize=25)
plt.ylabel(r'$\tan\beta$',fontsize=25)
plt.yscale('log')
if type == 1:
  plt.xlim([-0.601,0.601])
  plt.text(0.11, 0.2, r'2HDM Type-I', fontsize=15)
  plt.text(0.08, 0.16, r'(Run 2, 36/fb, ATLAS+CMS)', fontsize=12)
if type == 2:
  plt.xlim([-0.151,0.601])
  plt.text(0.21, 0.2, r'2HDM Type-II', fontsize=15)
  plt.text(0.21, 0.16, r'(Run 2, 36/fb, ATLAS+CMS)', fontsize=12)

fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "/results")
print("***** done *****")


