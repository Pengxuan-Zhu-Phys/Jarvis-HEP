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
from math import sqrt
from scipy.interpolate import UnivariateSpline

wdir = '/'.join(os.path.realpath(__file__).split("/")[:-1])+'/Grids/'

def gg_decay_ff(spline_deg=3):
    """read h -> g g partial widths grid @ BEST-QCD and interpolate"""
    
    GGfile = open(wdir+'GG_grid.dat',"r")
    GG_grid = {"TT": [], "CC": [], "TC": [], "BB": [], "TB": [], "CB": []}
    hmass = []

    for line in GGfile:
        line = line.strip("\n").split()
        hmass.append(float(line[0]))
        GG_grid["TT"].append(float(line[1]))
        GG_grid["CC"].append(float(line[2]))
        GG_grid["TC"].append(float(line[3]))
        GG_grid["BB"].append(float(line[4]))
        GG_grid["TB"].append(float(line[5]))
        GG_grid["CB"].append(float(line[6]))

    CggTT = UnivariateSpline(hmass, GG_grid["TT"], k=spline_deg, s=0)
    CggCC = UnivariateSpline(hmass, GG_grid["CC"], k=spline_deg, s=0)
    CggTC = UnivariateSpline(hmass, GG_grid["TC"], k=spline_deg, s=0)
    CggBB = UnivariateSpline(hmass, GG_grid["BB"], k=spline_deg, s=0)
    CggTB = UnivariateSpline(hmass, GG_grid["TB"], k=spline_deg, s=0)
    CggCB = UnivariateSpline(hmass, GG_grid["CB"], k=spline_deg, s=0)
    GGfile.close()

    gg_BESTQCD = {"CggTT":CggTT, "CggCC":CggCC, "CggTC":CggTC, "CggBB":CggBB,
                  "CggTB":CggTB, "CggCB":CggCB }

    return gg_BESTQCD



def gammagamma_ff(spline_deg=3):
    """read h -> gamma gamma partial widths grid @ BEST-QCD and interpolate """

    GaGafile = open(wdir+'GaGa_grid.dat',"r")
    gaga_grid = {"TT": [], "CC": [], "TC": [], "BB": [], "WW": [], "LL": [], "TB": [],
                 "CB": [], "TW": [], "CW": [], "BW": [], "TL": [], "CL": [], "BL": [], "LW": []}
    hmass = []

    for line in GaGafile:
        line = line.strip("\n").split()
        hmass.append(float(line[0]))
        gaga_grid["TT"].append(float(line[1]))
        gaga_grid["CC"].append(float(line[2]))
        gaga_grid["TC"].append(float(line[3]))
        gaga_grid["BB"].append(float(line[4]))
        gaga_grid["WW"].append(float(line[5]))
        gaga_grid["LL"].append(float(line[6]))
        gaga_grid["TB"].append(float(line[7]))
        gaga_grid["CB"].append(float(line[8]))
        gaga_grid["TW"].append(float(line[9]))
        gaga_grid["CW"].append(float(line[10]))
        gaga_grid["BW"].append(float(line[11]))
        gaga_grid["TL"].append(float(line[12]))
        gaga_grid["CL"].append(float(line[13]))
        gaga_grid["BL"].append(float(line[14]))
        gaga_grid["LW"].append(float(line[15]))

    CgagaTT = UnivariateSpline(hmass, gaga_grid["TT"], k=spline_deg, s=0)
    CgagaCC = UnivariateSpline(hmass, gaga_grid["CC"], k=spline_deg, s=0)
    CgagaTC = UnivariateSpline(hmass, gaga_grid["TC"], k=spline_deg, s=0)
    CgagaBB = UnivariateSpline(hmass, gaga_grid["BB"], k=spline_deg, s=0)
    CgagaWW = UnivariateSpline(hmass, gaga_grid["WW"], k=spline_deg, s=0)
    CgagaLL = UnivariateSpline(hmass, gaga_grid["LL"], k=spline_deg, s=0)
    CgagaTB = UnivariateSpline(hmass, gaga_grid["TB"], k=spline_deg, s=0)
    CgagaCB = UnivariateSpline(hmass, gaga_grid["CB"], k=spline_deg, s=0)
    CgagaTW = UnivariateSpline(hmass, gaga_grid["TW"], k=spline_deg, s=0)
    CgagaCW = UnivariateSpline(hmass, gaga_grid["CW"], k=spline_deg, s=0)
    CgagaBW = UnivariateSpline(hmass, gaga_grid["BW"], k=spline_deg, s=0)
    CgagaTL = UnivariateSpline(hmass, gaga_grid["TL"], k=spline_deg, s=0)
    CgagaCL = UnivariateSpline(hmass, gaga_grid["CL"], k=spline_deg, s=0)
    CgagaBL = UnivariateSpline(hmass, gaga_grid["BL"], k=spline_deg, s=0)
    CgagaLW = UnivariateSpline(hmass, gaga_grid["LW"], k=spline_deg, s=0)
    GaGafile.close()

    GaGa_BESTQCD = {"CgagaTT": CgagaTT ,"CgagaCC": CgagaCC ,"CgagaTC": CgagaTC ,"CgagaBB": CgagaBB ,"CgagaWW": CgagaWW ,"CgagaLL": CgagaLL ,"CgagaTB": CgagaTB ,"CgagaCB": CgagaCB ,"CgagaTW": CgagaTW ,"CgagaCW": CgagaCW ,"CgagaBW":CgagaBW ,"CgagaTL": CgagaTL ,"CgagaCL": CgagaCL ,"CgagaBL": CgagaBL ,"CgagaLW": CgagaLW}

    return GaGa_BESTQCD



def Zgamma_ff(spline_deg=3):
    """ read h -> Z gamma partial widths grid @ BEST-QCD and interpolate """

    ZGafile = open(wdir+'ZGa_grid.dat',"r")
    Zga_grid = {"TT": [], "CC": [], "TC": [], "BB": [], "WW": [], "LL": [], "TB": [],
                "CB": [], "TW": [], "CW": [], "BW": [], "TL": [], "CL": [], "BL": [], "LW": []}
    hmass = []

    for line in ZGafile:
        line = line.strip("\n").split()
        hmass.append(float(line[0]))
        Zga_grid["TT"].append(float(line[1]))
        Zga_grid["CC"].append(float(line[2]))
        Zga_grid["TC"].append(float(line[3]))
        Zga_grid["BB"].append(float(line[4]))
        Zga_grid["WW"].append(float(line[5]))
        Zga_grid["LL"].append(float(line[6]))
        Zga_grid["TB"].append(float(line[7]))
        Zga_grid["CB"].append(float(line[8]))
        Zga_grid["TW"].append(float(line[9]))
        Zga_grid["CW"].append(float(line[10]))
        Zga_grid["BW"].append(float(line[11]))
        Zga_grid["TL"].append(float(line[12]))
        Zga_grid["CL"].append(float(line[13]))
        Zga_grid["BL"].append(float(line[14]))
        Zga_grid["LW"].append(float(line[15]))

    CZgaTT = UnivariateSpline(hmass, Zga_grid["TT"], k=spline_deg, s=0)
    CZgaCC = UnivariateSpline(hmass, Zga_grid["CC"], k=spline_deg, s=0)
    CZgaTC = UnivariateSpline(hmass, Zga_grid["TC"], k=spline_deg, s=0)
    CZgaBB = UnivariateSpline(hmass, Zga_grid["BB"], k=spline_deg, s=0)
    CZgaWW = UnivariateSpline(hmass, Zga_grid["WW"], k=spline_deg, s=0)
    CZgaLL = UnivariateSpline(hmass, Zga_grid["LL"], k=spline_deg, s=0)
    CZgaTB = UnivariateSpline(hmass, Zga_grid["TB"], k=spline_deg, s=0)
    CZgaCB = UnivariateSpline(hmass, Zga_grid["CB"], k=spline_deg, s=0)
    CZgaTW = UnivariateSpline(hmass, Zga_grid["TW"], k=spline_deg, s=0)
    CZgaCW = UnivariateSpline(hmass, Zga_grid["CW"], k=spline_deg, s=0)
    CZgaBW = UnivariateSpline(hmass, Zga_grid["BW"], k=spline_deg, s=0)
    CZgaTL = UnivariateSpline(hmass, Zga_grid["TL"], k=spline_deg, s=0)
    CZgaCL = UnivariateSpline(hmass, Zga_grid["CL"], k=spline_deg, s=0)
    CZgaBL = UnivariateSpline(hmass, Zga_grid["BL"], k=spline_deg, s=0)
    CZgaLW = UnivariateSpline(hmass, Zga_grid["LW"], k=spline_deg, s=0)
    ZGafile.close()

    ZGa_BESTQCD = {"CZgaTT": CZgaTT ,"CZgaCC": CZgaCC ,"CZgaTC": CZgaTC ,"CZgaBB": CZgaBB ,"CZgaWW": CZgaWW ,"CZgaLL": CZgaLL ,"CZgaTB": CZgaTB ,"CZgaCB": CZgaCB ,"CZgaTW": CZgaTW ,"CZgaCW": CZgaCW ,"CZgaBW":CZgaBW ,"CZgaTL": CZgaTL ,"CZgaCL": CZgaCL ,"CZgaBL": CZgaBL ,"CZgaLW": CZgaLW}

    return ZGa_BESTQCD



def VBF_ff(spline_deg=3):
    """ read VBF -> h cross section (8 TeV) grid @ NLO-QCD and interpolate """

    VBF_NLO_file = open(wdir+'VBF_NLO8_grid.dat',"r")
    VBF_grid = {"WW": [], "ZZ": [], "WZ": []}
    hmassVBF = []

    for line in VBF_NLO_file:
        line = line.strip("\n").split()
        hmassVBF.append(float(line[0]))
        VBF_grid["WW"].append(float(line[1]))
        VBF_grid["ZZ"].append(float(line[2]))
        VBF_grid["WZ"].append(float(line[3])-float(line[1])-float(line[2]))

    CVBFW_NLO = UnivariateSpline(hmassVBF, VBF_grid["WW"], k=spline_deg)
    CVBFZ_NLO = UnivariateSpline(hmassVBF, VBF_grid["ZZ"], k=spline_deg)
    CVBFWZ_NLO = UnivariateSpline(hmassVBF, VBF_grid["WZ"], k=spline_deg)
    VBF_NLO_file.close()

    VBF_BESTQCD = {"CVBFW_NLO": CVBFW_NLO, "CVBFZ_NLO": CVBFZ_NLO, "CVBFWZ_NLO": CVBFWZ_NLO}

    return VBF_BESTQCD


def VBF13_ff(spline_deg=3):
    """ read VBF -> h cross section (13 TeV) grid @ NLO-QCD and interpolate """

    VBF13_NLO_file = open(wdir+'VBF_NLO13_grid.dat',"r")
    VBF13_grid = {"WW": [], "ZZ": [], "WZ": []}
    hmassVBF13 = []

    for line in VBF13_NLO_file:
        line = line.strip("\n").split()
        hmassVBF13.append(float(line[0]))
        VBF13_grid["WW"].append(float(line[1]))
        VBF13_grid["ZZ"].append(float(line[2]))
        VBF13_grid["WZ"].append(float(line[3])-float(line[1])-float(line[2]))

    CVBF13W_NLO = UnivariateSpline(hmassVBF13, VBF13_grid["WW"], k=spline_deg)
    CVBF13Z_NLO = UnivariateSpline(hmassVBF13, VBF13_grid["ZZ"], k=spline_deg)
    CVBF13WZ_NLO = UnivariateSpline(hmassVBF13, VBF13_grid["WZ"], k=spline_deg)
    VBF13_NLO_file.close()

    VBF13_BESTQCD = {"CVBF13W_NLO": CVBF13W_NLO, "CVBF13Z_NLO": CVBF13Z_NLO, "CVBF13WZ_NLO": CVBF13WZ_NLO}

    return VBF13_BESTQCD



def gg_prod_lhc8_ff(spline_deg=3):
    """ read g g -> h cross section (8TeV) grid @ NLO-QCD @ LHC8 and interpolate """

    ggF_NNLO_LHC8_file = open(wdir+'ggF_NNLO_LHC8_grid.dat',"r")
    ggF_LHC_grid = {"TT": [], "BB": [], "TB": []}
    hmassggF = []

    for line in ggF_NNLO_LHC8_file:
        line = line.strip("\n").split()
        hmassggF.append(float(line[0]))
        ggF_LHC_grid["TT"].append(float(line[1]))
        ggF_LHC_grid["BB"].append(float(line[2]))
        ggF_LHC_grid["TB"].append(float(line[3]))

    CggFT_NNLO_LHC8 = UnivariateSpline(hmassggF, ggF_LHC_grid["TT"], k=spline_deg, s=0)
    CggFB_NNLO_LHC8 = UnivariateSpline(hmassggF, ggF_LHC_grid["BB"], k=spline_deg, s=0)
    CggFTB_NNLO_LHC8 = UnivariateSpline(hmassggF, ggF_LHC_grid["TB"], k=spline_deg, s=0)
    ggF_NNLO_LHC8_file.close()

    ggF_LHC_BESTQCD = {"CggFT_NNLO_LHC8": CggFT_NNLO_LHC8, "CggFB_NNLO_LHC8": CggFB_NNLO_LHC8, "CggFTB_NNLO_LHC8": CggFTB_NNLO_LHC8}

    return ggF_LHC_BESTQCD


def gg_prod_lhc13_ff(spline_deg=3):
    """ read g g -> h cross section (13TeV) grid @ NLO-QCD @ LHC13 and interpolate """

    ggF_NNLO_LHC13_file = open(wdir+'ggF_NNLO_LHC13_grid.dat',"r")
    ggF_LHC_grid = {"TT": [], "BB": [], "TB": []}
    hmassggF = []

    for line in ggF_NNLO_LHC13_file:
        line = line.strip("\n").split()
        hmassggF.append(float(line[0]))
        ggF_LHC_grid["TT"].append(float(line[1]))
        ggF_LHC_grid["BB"].append(float(line[2]))
        ggF_LHC_grid["TB"].append(float(line[3]))

    CggFT_NNLO_LHC13 = UnivariateSpline(hmassggF, ggF_LHC_grid["TT"], k=spline_deg, s=0)
    CggFB_NNLO_LHC13 = UnivariateSpline(hmassggF, ggF_LHC_grid["BB"], k=spline_deg, s=0)
    CggFTB_NNLO_LHC13 = UnivariateSpline(hmassggF, ggF_LHC_grid["TB"], k=spline_deg, s=0)
    ggF_NNLO_LHC13_file.close()

    ggF13_LHC_BESTQCD = {"CggFT_NNLO_LHC13": CggFT_NNLO_LHC13, "CggFB_NNLO_LHC13": CggFB_NNLO_LHC13, "CggFTB_NNLO_LHC13": CggFTB_NNLO_LHC13}

    return ggF13_LHC_BESTQCD


def ggF_Tev_ff(spline_deg=3):
    """ read g g -> h cross section grid @ NLO-QCD @ Tevatron and interpolate """

    ggF_NNLO_Tev_file = open(wdir+'ggF_NNLO_Tev_grid.dat',"r")
    ggF_Tev_grid = {"TT": [], "BB": [], "TB": []}
    hmassggF = []

    for line in ggF_NNLO_Tev_file:
        line = line.strip("\n").split()
        hmassggF.append(float(line[0]))
        ggF_Tev_grid["TT"].append(float(line[1]))
        ggF_Tev_grid["BB"].append(float(line[2]))
        ggF_Tev_grid["TB"].append(float(line[3]))

    CggFT_NNLO_Tev = UnivariateSpline(hmassggF, ggF_Tev_grid["TT"], k=spline_deg, s=0)
    CggFB_NNLO_Tev = UnivariateSpline(hmassggF, ggF_Tev_grid["BB"], k=spline_deg, s=0)
    CggFTB_NNLO_Tev = UnivariateSpline(hmassggF, ggF_Tev_grid["TB"], k=spline_deg, s=0)
    ggF_NNLO_Tev_file.close()

    ggF_Tev_BESTQCD = {"CggFT_NNLO_Tev": CggFT_NNLO_Tev, "CggFB_NNLO_Tev": CggFB_NNLO_Tev, "CggFTB_NNLO_Tev": CggFTB_NNLO_Tev}

    return ggF_Tev_BESTQCD



#### decay: h -> gamma gamma @ BEST-QCD & reduced coupling ####

def redCgammagamma(CT, CB, CC, CL, CW, grid_interp, Cgammagammaadd=0.):
    gagaTT = grid_interp["CgagaTT"]
    gagaCC = grid_interp["CgagaCC"]
    gagaBB = grid_interp["CgagaBB"]
    gagaLL = grid_interp["CgagaLL"]
    gagaWW = grid_interp["CgagaWW"]
    gagaTB = grid_interp["CgagaTB"]
    gagaCB = grid_interp["CgagaCB"]
    gagaTL = grid_interp["CgagaTL"]
    gagaCL = grid_interp["CgagaCL"]
    gagaTW = grid_interp["CgagaTW"]
    gagaCW = grid_interp["CgagaCW"]
    gagaBW = grid_interp["CgagaBW"]
    gagaLW = grid_interp["CgagaLW"]
    gagaBL = grid_interp["CgagaBL"]
    gagaTC = grid_interp["CgagaTC"]
    
    amp_gaga_new = max(0.,(CT**2*gagaTT + CC**2*gagaCC + CB**2*gagaBB +
                    CL**2*gagaLL + CW**2*gagaWW +
                    CT*CB*gagaTB + CC*CB*gagaCB +
                    CT*CL*gagaTL + CC*CL*gagaCL +
                    CT*CW*gagaTW + CC*CW*gagaCW +
                    CB*CW*gagaBW + CL*CW*gagaLW +
                    CB*CL*gagaBL + CT*CC*gagaTC))
                                        
    amp_gaga = (gagaTT + gagaCC + gagaBB + gagaLL + gagaWW +
                gagaTB + gagaCB + gagaTL + gagaCL +
                gagaTW + gagaCW + gagaBW + gagaLW +
                gagaBL + gagaTC)
    
    return sqrt(amp_gaga_new/amp_gaga) + Cgammagammaadd


#### decay: h -> Z gamma @ BEST-QCD & reduced coupling ####

def redCZgamma(CT, CB, CC, CL, CW, grid_interp, CZgammaadd=0.):
    ZgaTT = grid_interp["CZgaTT"]
    ZgaCC = grid_interp["CZgaCC"]
    ZgaBB = grid_interp["CZgaBB"]
    ZgaLL = grid_interp["CZgaLL"]
    ZgaWW = grid_interp["CZgaWW"]
    ZgaTB = grid_interp["CZgaTB"]
    ZgaCB = grid_interp["CZgaCB"]
    ZgaTL = grid_interp["CZgaTL"]
    ZgaCL = grid_interp["CZgaCL"]
    ZgaTW = grid_interp["CZgaTW"]
    ZgaCW = grid_interp["CZgaCW"]
    ZgaBW = grid_interp["CZgaBW"]
    ZgaLW = grid_interp["CZgaLW"]
    ZgaBL = grid_interp["CZgaBL"]
    ZgaTC = grid_interp["CZgaTC"]
    
    amp_Zga_new = max(0.,(CT**2*ZgaTT + CC**2*ZgaCC + CB**2*ZgaBB +
                    CL**2*ZgaLL + CW**2*ZgaWW +
                    CT*CB*ZgaTB + CC*CB*ZgaCB +
                    CT*CL*ZgaTL + CC*CL*ZgaCL +
                    CT*CW*ZgaTW + CC*CW*ZgaCW +
                    CB*CW*ZgaBW + CL*CW*ZgaLW +
                    CB*CL*ZgaBL + CT*CC*ZgaTC) )
                    
    amp_Zga = (ZgaTT + ZgaCC + ZgaBB + ZgaLL + ZgaWW +
                ZgaTB + ZgaCB + ZgaTL + ZgaCL +
                ZgaTW + ZgaCW + ZgaBW + ZgaLW +
                ZgaBL + ZgaTC)
                    
    return sqrt(amp_Zga_new/amp_Zga)+CZgammaadd


#### decay: h -> g g @ BEST-QCD & reduced coupling ####

def redCgg(CT, CB, CC, grid_interp, Cggadd=0.):
    ggTT = grid_interp["CggTT"]
    ggCC = grid_interp["CggCC"]
    ggBB = grid_interp["CggBB"]
    ggTB = grid_interp["CggTB"]
    ggCB = grid_interp["CggCB"]
    ggTC = grid_interp["CggTC"]
    
    return sqrt( (CT**2*ggTT + CC**2*ggCC + CB**2*ggBB +
            CT*CB*ggTB + CC*CB*ggCB + CT*CC*ggTC)/
          (ggTT + ggCC + ggBB + ggTB + ggCB + ggTC) )


#### 8 TeV production: VBF -> h @ NLO-QCD & reduced coupling ####

def redCVBF(CW, CZ, grid_interp):
    VBFW_NLO =  grid_interp["CVBFW_NLO"]
    VBFZ_NLO =  grid_interp["CVBFZ_NLO"]
    VBFWZ_NLO = grid_interp["CVBFWZ_NLO"]
    
    return sqrt( (CW**2*VBFW_NLO + CZ**2*VBFZ_NLO + CW*CZ*VBFWZ_NLO)/
                  (VBFW_NLO + VBFZ_NLO + VBFWZ_NLO) )


#### 13 TeV production: VBF -> h @ NLO-QCD & reduced coupling ####

def redCVBF13(CW, CZ, grid_interp):
    VBFW_NLO =  grid_interp["CVBF13W_NLO"]
    VBFZ_NLO =  grid_interp["CVBF13Z_NLO"]
    VBFWZ_NLO = grid_interp["CVBF13WZ_NLO"]

    return sqrt( (CW**2*VBFW_NLO + CZ**2*VBFZ_NLO + CW*CZ*VBFWZ_NLO)/
                  (VBFW_NLO + VBFZ_NLO + VBFWZ_NLO) )


#### production: g g -> h @ BEST-QCD & reduced coupling ####

def redCggF_LHC8(CT, CB, grid_interp):
    ggFT_NNLO_LHC8 =  grid_interp["CggFT_NNLO_LHC8"]
    ggFB_NNLO_LHC8 =  grid_interp["CggFB_NNLO_LHC8"]
    ggFTB_NNLO_LHC8 = grid_interp["CggFTB_NNLO_LHC8"]
    
    return sqrt( (CT**2*ggFT_NNLO_LHC8 + CB**2*ggFB_NNLO_LHC8 + CT*CB*ggFTB_NNLO_LHC8)/
                 (ggFT_NNLO_LHC8 + ggFB_NNLO_LHC8 + ggFTB_NNLO_LHC8) )


#### production: g g -> h @ BEST-QCD @ LHC13 & reduced coupling ####

def redCggF_LHC13(CT, CB, grid_interp):
    ggFT_NNLO_LHC13 =  grid_interp["CggFT_NNLO_LHC13"]
    ggFB_NNLO_LHC13 =  grid_interp["CggFB_NNLO_LHC13"]
    ggFTB_NNLO_LHC13 = grid_interp["CggFTB_NNLO_LHC13"]

    return sqrt( (CT**2*ggFT_NNLO_LHC13 + CB**2*ggFB_NNLO_LHC13 + CT*CB*ggFTB_NNLO_LHC13)/
                 (ggFT_NNLO_LHC13 + ggFB_NNLO_LHC13 + ggFTB_NNLO_LHC13) )


#### production: g g -> h @ BEST-QCD @ Tevatron & reduced coupling ####

def redCggF_Tev(CT, CB, grid_interp):
    ggFT_NNLO_Tev =  grid_interp["CggFT_NNLO_Tev"]
    ggFB_NNLO_Tev =  grid_interp["CggFB_NNLO_Tev"]
    ggFTB_NNLO_Tev = grid_interp["CggFTB_NNLO_Tev"]
    
    
    return sqrt( (CT**2*ggFT_NNLO_Tev + CB**2*ggFB_NNLO_Tev + CT*CB*ggFTB_NNLO_Tev)/
                 (ggFT_NNLO_Tev + ggFB_NNLO_Tev + ggFTB_NNLO_Tev) )

# begin LDN added
# ref: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG2KAPPA
#### 8 TeV production, mH = 125 GeV: pp -> tHq (t-channel) cross section & reduced coupling ####
def redCtHq(CW, CT):
    return sqrt( 2.984*CT**2 + 3.886*CW**2 - 5.870*CT*CW )

#### 13 TeV production, mH = 125 GeV: pp -> tHq (t-channel) cross section & reduced coupling ####
def redCtHq13(CW, CT):
    return sqrt( 2.633*CT**2 + 3.578*CW**2 - 5.211*CT*CW )

#### 8 TeV production, mH = 125 GeV: pp -> tHW cross section & reduced coupling ####
def redCtHW(CW, CT):
    return sqrt( 2.426*CT**2 + 1.818*CW**2 - 3.244*CT*CW )

#### 13 TeV production, mH = 125 GeV: pp -> tHW cross section & reduced coupling ####
def redCtHW13(CW, CT):
    return sqrt( 2.909*CT**2 + 2.310*CW**2 - 4.220*CT*CW )

#### 8 TeV production, mH = 125 GeV: gg -> ZH cross section & reduced coupling ####
def redCggZH(CZ, CT, CB):
    return sqrt( 0.372*CT**2 + 0.0004*CB**2 + 2.302*CZ**2 + 0.003*CT*CB - 1.663*CT*CZ - 0.013*CB*CZ )

#### 13 TeV production, mH = 125 GeV: gg -> ZH cross section & reduced coupling ####
def redCggZH13(CZ, CT, CB):
    return sqrt( 0.456*CT**2 + 0.0004*CB**2 + 2.455*CZ**2 + 0.003*CT*CB - 1.902*CT*CZ - 0.011*CB*CZ )
# end LDN added
