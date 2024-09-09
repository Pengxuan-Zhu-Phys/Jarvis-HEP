##########################################################################
#
#  This file is part of Lilith
#  v1 (2015) by Jeremy Bernon and Beranger Dumont 
#  v2 (2019) by Sabine Kraml, Tran Quang Loc, Dao Thi Nhung, Le Duc Ninh 
#            converted to Python 3 by Marius Bertrand (Jul/Aug 2020)
#
#  Web page: http://lpsc.in2p3.fr/projects-th/lilith/
#
#  In case of questions email sabine.kraml@lpsc.in2p3.fr 
#
#
#    Lilith is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Lilith is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Lilith.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

import os
import numpy as np
from math import sqrt, log
from cmath import sqrt as csqrt
from cmath import asin as casin
from cmath import log as clog
from scipy.interpolate import UnivariateSpline
from .param import *

wdir = '/'.join(os.path.realpath(__file__).split("/")[:-1])+'/Grids/'

#### read VBF -> h cross section grid @ LO and interpolation ####

def VBF_ff(spline_deg=3):

    VBF_LO_file = open(wdir+'VBF_LO8_grid.dat',"r")
    VBF_LO_grid = {"WW": [], "ZZ": [], "WZ": []}
    hmassVBF = []
    for line in VBF_LO_file:
        line = line.strip("\n").split()
        hmassVBF.append(float(line[0]))
        VBF_LO_grid["WW"].append(float(line[1]))
        VBF_LO_grid["ZZ"].append(float(line[2]))
        VBF_LO_grid["WZ"].append(float(line[3])-float(line[1])-float(line[2]))

    CVBFW_LO = UnivariateSpline(hmassVBF, VBF_LO_grid["WW"], k=spline_deg, s=0)
    CVBFZ_LO = UnivariateSpline(hmassVBF, VBF_LO_grid["ZZ"], k=spline_deg, s=0)
    CVBFWZ_LO = UnivariateSpline(hmassVBF, VBF_LO_grid["WZ"], k=spline_deg, s=0)
    VBF_LO_file.close()

    VBF_LO = {"CVBFW_LO": CVBFW_LO, "CVBFZ_LO": CVBFZ_LO, "CVBFWZ_LO": CVBFWZ_LO}

    return VBF_LO


def VBF13_ff(spline_deg=3):

    VBF13_LO_file = open(wdir+'VBF_LO13_grid.dat',"r")
    VBF13_LO_grid = {"WW": [], "ZZ": [], "WZ": []}
    hmassVBF13 = []
    for line in VBF13_LO_file:
        line = line.strip("\n").split()
        hmassVBF13.append(float(line[0]))
        VBF13_LO_grid["WW"].append(float(line[1]))
        VBF13_LO_grid["ZZ"].append(float(line[2]))
        VBF13_LO_grid["WZ"].append(float(line[3])-float(line[1])-float(line[2]))

    CVBF13W_LO = UnivariateSpline(hmassVBF13, VBF13_LO_grid["WW"], k=spline_deg, s=0)
    CVBF13Z_LO = UnivariateSpline(hmassVBF13, VBF13_LO_grid["ZZ"], k=spline_deg, s=0)
    CVBF13WZ_LO = UnivariateSpline(hmassVBF13, VBF13_LO_grid["WZ"], k=spline_deg, s=0)
    VBF13_LO_file.close()

    VBF13_LO = {"CVBF13W_LO": CVBF13W_LO, "CVBF13Z_LO": CVBF13Z_LO, "CVBF13WZ_LO": CVBF13WZ_LO}

    return VBF13_LO


def fhiggs(t):
    if t<=1.:
        return casin(sqrt(t))**2.
    else:
        return -(log((sqrt(t) + sqrt(t-1.))/(sqrt(t) - sqrt(t-1.))) - pi*1j )**2./4.

def ghiggs(t):
    if t<=1:
        return csqrt(1-1/t)/2. * ( clog((1 + csqrt(1-1/t))/(1 - csqrt(1-1/t)))-pi*1j)
    else:
        return csqrt(1/t-1)*casin(csqrt(t))


def I1(tau,l):
    return (tau*l/(2.*(tau-l)) + tau**2*l**2 /(2.*(tau-l)**2)*(fhiggs(1/tau)-fhiggs(1/l)) +
           tau**2 * l /(tau-l)**2 * (ghiggs(1/tau)-ghiggs(1/l)))

def I2(tau,l):
    return -tau*l/(2.*(tau-l))*(fhiggs(1/tau)-fhiggs(1/l))


def A12(tau):
    return 2./tau *(1.+(1.-1./tau) * fhiggs(tau))

def A1(tau):
    return -(3.*tau+2.*tau**2. +3.*(2.*tau-1.) * fhiggs(tau))/tau**2

def A12Zgamma(tau,l):
    return I1(tau,l)-I2(tau,l)

def A1Zgamma(tau,l):
    return cW*(4.*(3.-sW2/cW2)*I2(tau,l)+((1.+2./tau)*sW2/cW2-(5.+2./tau))*I1(tau,l))

def A12A(tau):
    return 2/tau*fhiggs(tau)


def computeformfactors():
    FF = {}
    
    FF["A12t"] = lambda mh: A12((mh/(2.*mt))**2)
    FF["A12c"] = lambda mh: A12((mh/(2.*mc))**2)
    FF["A12b"] = lambda mh: A12((mh/(2.*mb))**2)
    FF["A12tau"] = lambda mh: A12((mh/(2.*mtau))**2)
    FF["A1W"] = lambda mh: A1((mh/(2.*mW))**2)
    FF["A12At"] = lambda mh: A12A((mh/(2.*mt))**2)
    FF["A12Ac"] = lambda mh: A12A((mh/(2.*mc))**2)
    FF["A12Ab"] = lambda mh: A12A((mh/(2.*mb))**2)
    FF["A12Atau"] = lambda mh: A12A((mh/(2.*mtau))**2)

    
    FF["A12Zt"] = lambda mh: A12Zgamma(4*(mt/(mh*1.))**2, 4*(mt/(mZ*1.))**2)
    FF["A12Zc"] = lambda mh: A12Zgamma(4*(mc/(mh*1.))**2, 4*(mc/(mZ*1.))**2)
    FF["A12Zb"] = lambda mh: A12Zgamma(4*(mb/(mh*1.))**2, 4*(mb/(mZ*1.))**2)
    FF["A12Ztau"] = lambda mh: A12Zgamma(4*(mtau/(mh*1.))**2, 4*(mtau/(mZ*1.))**2)
    FF["A1ZW"] = lambda mh: A1Zgamma(4*(mW/(mh*1.))**2, 4*(mW/(mZ*1.))**2)
    FF["A12AZt"] = lambda mh: I2(4*(mt/(mh*1.))**2, 4*(mt/(mZ*1.))**2)
    FF["A12AZc"] = lambda mh: I2(4*(mc/(mh*1.))**2, 4*(mc/(mZ*1.))**2)
    FF["A12AZb"] = lambda mh: I2(4*(mb/(mh*1.))**2, 4*(mb/(mZ*1.))**2)
    FF["A12AZtau"] = lambda mh: I2(4*(mtau/(mh*1.))**2, 4*(mtau/(mZ*1.))**2)

    return FF

#### decay: h -> gamma gamma width @ LO & reduced coupling ####
def Htogammagamma(mh, CT, CB, CC, CL, CW, CTIM, CBIM, CCIM, CLIM, FF):
    return      (10**6*Gf*alpha**2/(128.*pi**3*np.sqrt(2))*mh**3 *
                abs(3.*(2./3.)**2 *(CT*FF["A12t"] + CC*FF["A12c"]) +
                (CB*3.*(1./3.)**2 * FF["A12b"] + CL*FF["A12tau"])+CW*FF["A1W"])**2. +
                10**6*Gf*alpha**2/(128.*pi**3*np.sqrt(2))*mh**3 *
                abs(3.*(2./3.)**2 *(CTIM*FF["A12At"] + CCIM*FF["A12Ac"]) +
                (CBIM*3.*(1./3.)**2 * FF["A12Ab"] + CLIM*FF["A12Atau"]))**2.)


def redCgammagamma(CT, CB, CC, CL, CW, CTIM, CBIM, CCIM, CLIM, FF, Cgammagammaadd=0.):
    A12t = FF["A12t"]
    A12c = FF["A12c"]
    A12b = FF["A12b"]
    A12tau = FF["A12tau"]
    A1W = FF["A1W"]
    A12At = FF["A12At"]
    A12Ac = FF["A12Ac"]
    A12Ab = FF["A12Ab"]
    A12Atau = FF["A12Atau"]

    return (sqrt( ( (abs(3.*(2./3.)**2 *(CT*A12t + CC*A12c) +
                       CB*3.*(1./3.)**2 * A12b + CL*A12tau+CW*A1W)**2.) +
                   (abs(3.*(2./3.)**2 *(CTIM*A12At + CCIM*A12Ac) +
                        3.*(-1./3.)**2*CBIM*A12Ab + CLIM*A12Atau)**2)) /
                 (abs(3.*(2./3.)**2 *(A12t + A12c) +
                     (3.*(1./3.)**2 * A12b + A12tau)+A1W)**2.) )
                + Cgammagammaadd)


#### decay: h -> Z gamma width @ LO & reduced coupling ####
def HtoZgamma(mh, CT, CB, CC, CL, CW, CTIM, CBIM, CCIM, CLIM, FF):
    return (10**6*Gf**2*mW**2*alpha*mh**3/(64.*pi**4)*(1-mZ**2/mh**2)**3 *
        abs( 1/(cW)*3.*2/3.*(CT*(2*1/2. - 4*2/3.*sW2)*FF["A12Zt"] +
                             CC*(2*1/2. - 4*2/3.*sW2)*FF["A12Zc"]) +
            1/(cW)*(3*(-1/3.)*CB*(2*(-1/2.) - 4*(-1/3.)*sW2)*FF["A12Zb"] +
                    (-1)*CL*(2*(-1/2.) - 4*(-1)*sW2)*FF["A12Ztau"]) +
            CW*FF["A1ZW"] )**2 +
            10**6*Gf**2*mW**2*alpha*mh**3/(16.*pi**4)*(1-mZ**2/mh**2)**3 *
                abs( 1/(cW)*3.*2/3.*(CTIM*(2*1/2. - 4*2/3.*sW2)*FF["A12AZt"] +
                                     CCIM*(2*1/2. - 4*2/3.*sW2)*FF["A12AZc"]) +
                    1/(cW)*(3*(-1/3.)*CBIM*(2*(-1/2.) - 4*(-1/3.)*sW2)*FF["A12AZb"] +
                            (-1)*CL*(2*(-1/2.) - 4*(-1)*sW2)*FF["A12Ztau"]) )**2)


def redCZgamma(CT, CB, CC, CL, CW, CTIM, CBIM, CCIM, CLIM, FF, CZgammaadd=0.):
    A12Zt = FF["A12Zt"]
    A12Zc = FF["A12Zc"]
    A12Zb = FF["A12Zb"]
    A12Ztau = FF["A12Ztau"]
    A1ZW = FF["A1ZW"]
    A12AZt = FF["A12AZt"]
    A12AZc = FF["A12AZc"]
    A12AZb = FF["A12AZb"]
    A12AZtau = FF["A12AZtau"]
    vt = (2*1/2. - 4*2/3.*sW2)
    vc = (2*1/2. - 4*2/3.*sW2)
    vb = (2*(-1/2.) - 4*(-1/3.)*sW2)
    vl = (2*(-1/2.) - 4*(-1)*sW2)
    
    return (sqrt( (abs( 1/(cW)*(3.*2/3.*(CT*vt*A12Zt + CC*vc*A12Zc) +
                       (3*(-1/3.)*CB*vb*A12Zb +
                        (-1)*CL*vl*A12Ztau)) + CW*A1ZW )**2 +
                  4*abs(1/(cW)*(3.*2/3.*(CTIM*vt*A12AZt + CCIM*vc*A12AZc) +
                        3*(-1/3.)*CBIM*vb*A12AZb + (-1)*CLIM*vl*A12AZtau))**2)/
                (abs(1/(cW)*(3.*2/3.*(vt*A12Zt + vc*A12Zc) +
                            (3*(-1/3.)*vb*A12Zb + (-1)*vl*A12Ztau)) + A1ZW )**2) )
            + CZgammaadd)


#### decay: h -> g g width @ LO & reduced coupling ####
def Htogg(mh, CT, CB, CC, CTIM, CBIM, CCIM, FF):
    return (10**3*Gf*alphas**2*mh**3/(36.*np.sqrt(2)*pi**3) *
        abs(0.75*(CT*FF["A12t"] + CB*FF["A12b"] + CC*FF["A12c"]))**2 +
            10**3*Gf*alphas**2*mh**3/(36.*np.sqrt(2)*pi**3) *
            abs(0.75*(CTIM*FF["A12At"] + CBIM*FF["A12Ab"] + CCIM*FF["A12Ac"]))**2)


def redCgg(CT, CB, CC, CTIM, CBIM, CCIM, FF, Cggadd=0.):
    A12t = FF["A12t"]
    A12c = FF["A12c"]
    A12b = FF["A12b"]
    A12At = FF["A12At"]
    A12Ac = FF["A12Ac"]
    A12Ab = FF["A12Ab"]
    
    return (sqrt( (abs(0.75*(CT*A12t + CB*A12b + CC*A12c))**2 +
                  abs(0.75*(CTIM*A12At + CBIM*A12Ab + CCIM*A12Ac))**2)/
                 (abs(0.75*(A12t + A12b + A12c))**2) )
           + Cggadd)


#### production: g g -> h cross section @ LO ####
def ggFh(mh, CT, CB, CC, CTIM, CBIM, CCIM, FF):
    return (Gf*alphas_mh**2/(288.*np.sqrt(2)*pi) *
        abs(0.75*(CT*FF["A12t"] + CB*FF["A12b"] + CC*FF["A12c"]))**2 +
            Gf*alphas_mh**2/(288.*np.sqrt(2)*pi) *
        abs(0.75*(CTIM*FF["A12At"] + CBIM*FF["A12Ab"] + CCIM*FF["A12Ac"]))**2)


#### 8 TeV production: VBF -> h cross section @ LO & reduced coupling ####
def redCVBF(CW, CZ, grid_interp):
    VBFW_LO = grid_interp["CVBFW_LO"]
    VBFZ_LO = grid_interp["CVBFZ_LO"]
    VBFWZ_LO = grid_interp["CVBFWZ_LO"]
    
    return sqrt( (CW**2*VBFW_LO + CZ**2*VBFZ_LO + CW*CZ*VBFWZ_LO)/
                 (VBFW_LO + VBFZ_LO + VBFWZ_LO) )

#### 13 TeV production: VBF -> h cross section @ LO & reduced coupling ####
def redCVBF13(CW, CZ, grid_interp):
    VBFW_LO = grid_interp["CVBF13W_LO"]
    VBFZ_LO = grid_interp["CVBF13Z_LO"]
    VBFWZ_LO = grid_interp["CVBF13WZ_LO"]

    return sqrt( (CW**2*VBFW_LO + CZ**2*VBFZ_LO + CW*CZ*VBFWZ_LO)/
                 (VBFW_LO + VBFZ_LO + VBFWZ_LO) )
