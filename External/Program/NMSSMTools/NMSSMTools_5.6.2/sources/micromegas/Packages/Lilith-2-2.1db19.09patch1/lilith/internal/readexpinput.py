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

import sys
try:
    from lxml import etree
except:
    import xml.etree.ElementTree as etree
from ..errors import ExpInputError, ExpInputIOError
from scipy import interpolate
import numpy as np
from scipy.optimize import fsolve
import math
from . import brsm as BR_SM
from warnings import warn

class ReadExpInput:
    """Read the experimental input in XML and extracts all information."""

    def __init__(self):
        """Initialize the VBF, WH and ZH cross section ratios."""

        self.eff_VVH = None
        self.eff_top = None
        self.mu = []
        self.filepath = ""

    def warning(self, message):
        """Customized warnings."""

        warn("in the file " + self.filepath + ": " + message, Warning,
                      stacklevel=3)

    def get_filelist(self, filepath):
        """get list of files from .list experimental input"""

        expfiles = [] # list of all paths to XML input files

        filepath_split = filepath.split("/")
        expdata_dir = "/".join(filepath_split[:-1])

        try:
            with open(filepath) as linput:
                expfiles = [] # list of all paths to XML input files
                for line in linput:
                    # remove comment, space and new lines from the line
                    line = line.split("#")[0].rstrip("\n").strip()
                    if line == "": # empty line or comment
                        continue
                    if line[0] == "/": # absolute filepath
                        expfiles.append(line)
                    else: # relative filepath
                        expfiles.append(expdata_dir+"/"+line)
        except IOError as e:
            raise ExpInputIOError(
                'I/O error({0}): {1}'.format(e.errno, e.strerror) + '; cannot' +
                ' open the experimental list file "' + filepath + '".')
        
        return expfiles
    

    def read_file(self, filepath):
        """read individual xml files"""

        self.filepath = filepath

        root = self.produce_tree()
        if root.tag != "expmu":
            raise ExpInputError(self.filepath, "root tag is not <expmu>")

        (dim, decay, type) = self.get_mode(root)
        
        self.get_mass(root)
        (experiment, source, sqrts) = self.get_metadata(root)
       
        eff = self.read_eff(root, dim, decay)
        (bestfit, param, grid, Lxy, LChi2min) = self.read_mus(root, dim, type)

        self.mu.append({"filepath": self.filepath,
                        "dim": dim, "type": type,
                        "eff": eff,
                        "bestfit": bestfit, "param": param, "grid": grid,
                        "Lxy": Lxy, "LChi2min": LChi2min,
                        "experiment": experiment, "source": source,
                        "sqrts": sqrts, "eff": eff})

    def produce_tree(self):
        """Produce the XML tree with ElementTree."""

        try:
            with open(self.filepath) as f:
                tree = etree.parse(f)
        except IOError as e:
            raise ExpInputIOError(
                'I/O error({0}): {1}'.format(e.errno, e.strerror) + '; cannot' +
                'open the experimental input file "' + self.filepath + '".')

        return tree.getroot()

    def get_mode(self, root):
        """Get the dimension, decay and type of the experimental mu."""
        
        allowed_decays = ["gammagamma", "ZZ", "WW", "Zgamma",
                          "tautau", "bb", "cc", "mumu", "invisible", "gg"]

        mandatory_attribs = {"type":["n", "vn", "p", "f"]}  
                             # types vn and p added for variable Gaussian and Poisson fits

        optional_attribs = {"decay": allowed_decays}

        for mandatory_attrib, allowed_values in list(mandatory_attribs.items()):
            if mandatory_attrib not in root.attrib:
                # if "type" not in attribute
                raise ExpInputError(self.filepath,
                                    'mandatory attribute of root tag "' +
                                    mandatory_attrib + '" is not present.')
            
            if root.attrib[mandatory_attrib] not in allowed_values:
                # if type="z" for instance
                raise ExpInputError(self.filepath,
                                    'mandatory attribute of root tag "' +
                                    mandatory_attrib + '" has value "' +
                                   root.attrib[mandatory_attrib] +
                                   '" which is unknown. Allowed values are : ' +
                                   str(allowed_values))

        dim = int(root.attrib["dim"])
        if dim <= 0: 
            raise ExpInputError(self.filepath,
                                    'attribute dim is non-positive.')
        type = root.attrib["type"]

        decay = "mixture"

        for optional_attrib, allowed_values in list(optional_attribs.items()):
            if optional_attrib in root.attrib:
                # if "decay" in attribute
                if root.attrib[optional_attrib] not in allowed_values:
                    # if decay="yy" for instance
                    raise ExpInputError(self.filepath,
                                        'optional attribute of root tag "' +
                                        optional_attrib + '" has value "' +
                                        root.attrib[optional_attrib] +
                                        '" which is unknown. Allowed values ' +
                                        'are: ' + str(allowed_values))
                else:
                    decay = root.attrib["decay"]

        return (dim, decay, type)
    
    def get_mass(self, root):
        def_mass = 125. # default value
        mass = def_mass
    
        for child in root:
            if child.tag == "mass":
                try:
                    mass = float(child.text)
                except TypeError: # empty tag is of type NULL
                    self.warning('<mass> tag is empty; ' +
                                 'setting the mass to ' + def_mass + ' GeV')
                    mass = def_mass
                except ValueError:
                    raise ExpInputError(self.filepath,
                                        "value of <mass> tag is not a number")
        self.mass = mass

    def get_metadata(self, root):
        experiment = ""
        source = ""
        sqrts = ""

        for child in root:
            if child.tag == "experiment":
                experiment = child.text
            if child.tag == "source":
                source = child.text
            if child.tag == "sqrts":
                sqrts = child.text
            
        return (experiment, source, sqrts)

    def read_eff(self, root, dim, decay):
        allowed_decays = ["gammagamma", "ZZ", "WW", "Zgamma",
                          "tautau", "bb", "cc", "mumu", "invisible", "gg"]

        # read the efficiencies
        if dim == 1: # 1D signal strength
            eff = {"x": {}}
            axis_label = "x"

            mandatory_attribs = {"prod": ["ggH", "VVH", "VBF", "VH", "WH", "qqZH", "ggZH", "ZH", "ttH", "tHq", "tHW", "tH", "top", "bbH"]}
            if decay == "mixture":
                mandatory_attribs["decay"] = allowed_decays

            for child in root:
                if child.tag=="sqrts":
                    sqrts = child.text
                    if sqrts not in ["1.96","7","8","7.","8.","7.0","8.0","7+8","13","13.","13.0"]:
                        self.warning("sqrt(s) of experimental input is not a Tevatron (1.96) or 7,8,13 TeV LHC result."+" Lilith will use automatically the 13 TeV form factors and cross sections if needed.")
                if child.tag == "eff":
                    for mandatory_attrib, allowed_values in list(mandatory_attribs.items()):
                        if mandatory_attrib not in child.attrib:
                            # if "axis" or "prod" not in attribute
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" is not present.')
                        if child.attrib[mandatory_attrib] not in allowed_values:
                            # if axis="z" or prod="yy"
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" has value "' +
                                                child.attrib[mandatory_attrib] + '" which is unknown. Allowed values are : ' + str(allowed_values))

                    prod_label = child.attrib["prod"]
                    if decay == "mixture":
                        decay_label = child.attrib["decay"]
                    else:
                        decay_label = decay

                    if (prod_label,decay_label) in eff["x"]:
                        self.warning('<eff> tag with prod="' + prod_label +
                                     '" and decay="' + decay_label +
                                     '" is being redefined.')
                        
                    try:
                        eff[axis_label][prod_label,decay_label] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        self.warning('<eff> tag for axis="' + axis_label +
                                     '", prod="' + prod_label + '" and decay="' +
                                     decay_label + '" is empty; setting to ' +
                                     'default value of 0')
                        eff[axis_label][prod_label,decay_label] = 0.
                    except ValueError:
                        raise ExpInputError(self.filepath,
                                            'value of <eff> tag with axis="' + axis_label +
                                            '" and prod="' + prod_label + '" and decay="' + decay_label + '" is not a number')

        elif dim == 2: # 2D signal strength
            eff = {"x": {}, "y": {}}

            mandatory_attribs = {"axis": ["x", "y"],
                             "prod": ["ggH", "VVH", "VBF", "VH", "WH", "qqZH", "ggZH", "ZH", "ttH", "tHq", "tHW", "tH", "top", "bbH"]}
            if decay == "mixture":
                mandatory_attribs["decay"] = allowed_decays
            
            for child in root:
                if child.tag=="sqrts":
                        sqrts = child.text
                        if sqrts not in ["1.96","7","8","7.","8.","7.0","8.0","7+8","13","13.","13.0"]:
                            self.warning("sqrt(s) of experimental input is not a Tevatron (1.96) or 7,8,13 TeV LHC result."+
                            " Lilith will use automatically the 13 TeV form factors and cross sections if needed.")
                if child.tag == "eff":
                    for mandatory_attrib, allowed_values in list(mandatory_attribs.items()):
                        if mandatory_attrib not in child.attrib:
                            # if "axis" or "prod" not in attribute
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" is not present.')
                        if child.attrib[mandatory_attrib] not in allowed_values:
                            # if axis="z" or prod="yy"
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" has value "' +
                                                child.attrib[mandatory_attrib] + '" which is unknown. Allowed values are : ' + str(allowed_values))

                        axis_label = child.attrib["axis"]
                        prod_label = child.attrib["prod"]
                        if decay == "mixture":
                            decay_label = child.attrib["decay"]
                        else:
                            decay_label = decay


                    if (prod_label,decay_label) in eff[axis_label]:
                        self.warning('<eff> tag with axis="' + axis_label +
                                     '", prod="' + prod_label +
                                     '" and decay="' +decay_label +
                                     '" is being redefined.')
                        
                    try:
                        eff[axis_label][prod_label,decay_label] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        self.warning('<eff> tag for axis="' + axis_label +
                                     '", prod="' + prod_label + '" and decay="' +
                                     decay_label + '" is empty; setting to ' +
                                     'default value of 0')
                        eff[axis_label][prod_label,decay_label] = 0.
                    except ValueError:
                        raise ExpInputError(self.filepath,
                                            'value of <eff> tag with axis="' + axis_label +
                                            '" and prod="' + prod_label + '" and decay="' + decay_label + '" is not a number')
        ### new from v2.0 onwards 
        elif dim >= 3:
            eff = {"d1": {}, "d2": {}, "d3": {}}
            axis_attribs = ["d1", "d2", "d3"]
            for i in range(4,dim+1):
                d = "d" + str(i)
                eff.update({d: {}})
                axis_attribs.append(d)
            prod_attribs = ["ggH", "VVH", "VBF", "VH", "WH", "qqZH", "ggZH", "ZH", "ttH", "tHq", "tHW", "tH", "top", "bbH"]
            mandatory_attribs = {"axis": axis_attribs, "prod": prod_attribs}

            if decay == "mixture":
                mandatory_attribs["decay"] = allowed_decays

            for child in root:
                if child.tag=="sqrts":
                    sqrts = child.text
                    if sqrts not in ["1.96","7","8","7.","8.","7.0","8.0","7+8","13","13.","13.0"]:
                        self.warning("sqrt(s) of experimental input is not a Tevatron (1.96) or 7,8,13 TeV LHC result."+
                                " Lilith will use automatically the 13 TeV form factors and cross sections if needed.")
                if child.tag == "eff":
                    for mandatory_attrib, allowed_values in list(mandatory_attribs.items()):
                        if mandatory_attrib not in child.attrib:
                            # if "axis" or "prod" not in attribute
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" is not present.')
                        if child.attrib[mandatory_attrib] not in allowed_values:
                            # if axis="z" or prod="yy"
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" has value "' +
                                                child.attrib[mandatory_attrib] + '" which is unknown. Allowed values are : ' + str(allowed_values))

                        axis_label = child.attrib["axis"]
                        prod_label = child.attrib["prod"]
                        if decay == "mixture":
                            decay_label = child.attrib["decay"]
                        else:
                            decay_label = decay


                    if (prod_label,decay_label) in eff[axis_label]:
                        self.warning('<eff> tag with axis="' + axis_label +
                                     '", prod="' + prod_label +
                                     '" and decay="' +decay_label +
                                     '" is being redefined.')
                        
                    try:
                        eff[axis_label][prod_label,decay_label] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        self.warning('<eff> tag for axis="' + axis_label +
                                     '", prod="' + prod_label + '" and decay="' +
                                     decay_label + '" is empty; setting to ' +
                                     'default value of 0')
                        eff[axis_label][prod_label,decay_label] = 0.
                    except ValueError:
                        raise ExpInputError(self.filepath,
                                            'value of <eff> tag with axis="' + axis_label +
                                            '" and prod="' + prod_label + '" and decay="' + decay_label + '" is not a number')

	
        if sqrts in ["1.96","7","8","7.","8.","7.0","8.0","7+8"]:
            self.eff_VVH = BR_SM.geteffVVHfunctions(8)
            self.eff_top = BR_SM.getefftopfunctions(8)
        else:
            self.eff_VVH = BR_SM.geteffVVHfunctions(13)
            self.eff_top = BR_SM.getefftopfunctions(13)
		
        effWH_VH = self.eff_VVH["eff_WH"](self.mass) # relative to VH

        effqqZH_ZH = self.eff_VVH["eff_qqZH"](self.mass) # relative to ZH
        effggZH_ZH = self.eff_VVH["eff_ggZH"](self.mass) # relative to ZH
        
        effVBF_VVH = self.eff_VVH["eff_VBF"](self.mass) # relative to VVH

        effVH_VVH = self.eff_VVH["eff_VH"](self.mass) # relative to VVH
        effZH_VH = self.eff_VVH["eff_ZH"](self.mass) # relative to VH

        effqqZH_VH = effZH_VH*effqqZH_ZH
        effggZH_VH = effZH_VH*effggZH_ZH

        effWH_VVH = effVH_VVH*effWH_VH # relative to VVH
        effqqZH_VVH = effVH_VVH*effqqZH_VH # relative to VVH
        effggZH_VVH = effVH_VVH*effggZH_VH # relative to VVH
        effZH_VVH = effVH_VVH*effZH_VH # relative to VVH

        efftHq_tH = self.eff_top["eff_tHq"](self.mass) # relative to tH
        efftHW_tH = self.eff_top["eff_tHW"](self.mass) # relative to tH

        effttH_top = self.eff_top["eff_ttH"](self.mass) # relative to top
        efftH_top = self.eff_top["eff_tH"](self.mass) # relative to top
        efftHq_top = efftH_top*efftHq_tH # relative to top
        efftHW_top = efftH_top*efftHW_tH # relative to top

        multiprod = {"ZH": {"qqZH": effqqZH_ZH, "ggZH": effggZH_ZH}, "VH": {"WH": effWH_VH, "qqZH": effqqZH_VH, "ggZH": effggZH_VH}, "VVH": {"VBF": effVBF_VVH, "WH": effWH_VVH, "qqZH": effqqZH_VVH, "ggZH": effggZH_VVH}, "tH": {"tHq": efftHq_tH, "tHW": efftHW_tH}, "top": {"ttH": effttH_top, "tHq": efftHq_top, "tHW": efftHW_top}}
        
        if dim == 1:
            self.check_multiprod(eff["x"], multiprod)
        elif dim == 2:
            self.check_multiprod(eff["x"], multiprod)
            self.check_multiprod(eff["y"], multiprod)
        elif dim >= 3:
            for i in range(1,dim+1):
                self.check_multiprod(eff["d"+str(i)], multiprod)

        # now all reduced couplings have been properly defined, one can
        # delete all multiparticle labels
        if dim == 1:
            effCleanX = eff["x"].copy()
            for (p,decay) in eff["x"]:
                if p in multiprod:
                    del effCleanX[p,decay]
            eff["x"] = effCleanX
        elif dim == 2:
            effCleanX = eff["x"].copy()
            for (p,decay) in eff["x"]:
                if p in multiprod:
                    del effCleanX[p,decay]
            effCleanY = eff["y"].copy()
            for (p,decay) in eff["y"]:
                if p in multiprod:
                    del effCleanY[p,decay]
            eff["x"] = effCleanX
            eff["y"] = effCleanY
        elif dim >= 3:
            for i in range(1,dim+1):
                d = "d" + str(i)
                effCleanD = eff[d].copy()
                for (p,decay) in eff[d]:
                    if p in multiprod:
                        del effCleanD[p,decay]
                eff[d] = effCleanD

        # check that efficiencies add up to 1, otherwise issue a warning
        # or an error
        for axis in eff:
            sumeff = 0
            for prod in eff[axis]:
                sumeff += eff[axis][prod]
            
            if sumeff == 0:
                raise ExpInputError(self.filepath,
                                    "no <eff> tag found for " + axis + " axis")

            if sumeff < 0.99:
                self.warning('the sum of efficiencies for axis="' +
                             axis + '" is less than 1 (value: ' +
                             str(sumeff) + ')')
            elif sumeff > 1.01:
                raise ExpInputError(self.filepath,
                                    'the sum of efficiencies for axis="' +
                                    axis + '" is greater than 1 (value: ' +
                                    str(sumeff) + ')')

        return eff


    def read_mus(self, root, dim, type):
        # first, read the bestfit
        bestfit = {}
        LChi2min = 0
        
        for child in root:
            if child.tag == "bestfit":
                if type == "f":
                    self.warning('block <bestfit> in experimental mu of ' +
                                 'type "full"... skipping.')
                    continue

                if dim == 1:
                    # read directly the number
                    
                    if "x" in bestfit:
                        self.warning("redefinition of the bestfit...")
                    
                    try:
                        bestfit["x"] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        self.warning('<x> tag in <bestfit> block is empty; ' +
                                     'setting to 0')
                        bestfit["x"] = 0.
                    except ValueError:
                        raise ExpInputError(self.filepath,
                                            "value of <besfit> tag is not a number")
            
                elif dim == 2:
                    bestfit_allowedsubtags = ["x", "y"]
                
                    for bfit in child:
                        if bfit.tag in bestfit_allowedsubtags:
                            
                            if bfit.tag in bestfit:
                                self.warning("redefinition of the bestfit...")
                        
                            try:
                                bestfit[bfit.tag] = float(bfit.text)
                            except TypeError: # empty tag is of type NULL
                                self.warning('<' + bfit.tag + '> tag in ' +
                                             '<bestfit> block is empty; ' +
                                             'setting to 0')
                                bestfit[bfit.tag] = 0.
                            except ValueError:
                                raise ExpInputError(self.filepath,
                                                    "value of <besfit> tag is not a number")
                        else:
                            raise ExpInputError(self.filepath,
                                                "subtag in bestfit not known")
                elif dim >= 3:
                    bestfit_allowedsubtags = ["d1", "d2", "d3"]
                    for i in range(4,dim+1):
                        bestfit_allowedsubtags.append("d"+str(i))
                
                    for bfit in child:
                        if bfit.tag in bestfit_allowedsubtags:
                            
                            if bfit.tag in bestfit:
                                self.warning("redefinition of the bestfit...")
                        
                            try:
                                bestfit[bfit.tag] = float(bfit.text)
                            except TypeError: # empty tag is of type NULL
                                self.warning('<' + bfit.tag + '> tag in ' +
                                             '<bestfit> block is empty; ' +
                                             'setting to 0')
                                bestfit[bfit.tag] = 0.
                            except ValueError:
                                raise ExpInputError(self.filepath,
                                                    "value of <besfit> tag is not a number")
                        else:
                            raise ExpInputError(self.filepath,
                                                "subtag in bestfit not known")
                    
                if dim == 1 and "x" not in bestfit:
                    raise ExpInputError(self.filepath,
                                        "best fit point should be specified.")
                if dim == 2 and ("x" not in bestfit or "y" not in bestfit):
                    raise ExpInputError(self.filepath,
                                        "best fit point should be specified for x and y.")
                if dim >= 3:
                    for i in range(dim):
                        if (bestfit_allowedsubtags[i] not in bestfit):
                            raise ExpInputError(self.filepath,
                                        "best fit point should be specified for "+bestfit_allowedsubtags[i])
        
        # then, read the param...
        param = {}
        
        for child in root:
            if child == "param":
                break
        param_tag = child

        param["uncertainty"] = {}

### v2.0 onwards: 
# Poisson likelihood, type = "p", dim = 1 or 2 with corralation
# variable Gaussian, type = "vn", dim >= 1 with corralations

        if dim == 2:
            param["uncertainty"]["x"] = {}
            param["uncertainty"]["y"] = {}
        elif dim >= 3:
            axis_allowedtags = []
            for i in range(1,dim+1):
                d = "d" + str(i)
                axis_allowedtags.append(d)
                param["uncertainty"][d] = {}
            param["correlation"] = {}
            corr_allowedtags = []
            for i in range(1,dim+1):
                for j in range(i+1,dim+1):
                    dd = "d" + str(i) + "d" + str(j)
                    corr_allowedtags.append(dd)
        if type == "p" and dim == 2:
            param["gamma"] = {}
            param["nu"] = {}

        for child in param_tag:
            if child.tag is etree.Comment:
                # ignore all comments
                continue

            if dim == 1 and (type == "n" or type == "vn" or type == "p"):
                if child.tag == "uncertainty":
                    if "side" not in child.attrib:
                        try:
                            unc_value = float(child.text)
                        except TypeError: # empty tag is of type NULL
                            self.warning('<uncertainty> tag is empty; ' +
                                         'setting to 0')
                            unc_value = 0.
                        except ValueError:
                            raise ExpInputError(self.filepath,
                                                "value of <uncertainty> tag is not a number")

                        param["uncertainty"]["left"] = unc_value
                        param["uncertainty"]["right"] = unc_value
                    else:
                        if child.attrib["side"] not in ["left", "right"]:
                            raise ExpInputError(self.filepath,
                                                "attribute of uncertainty is not left nor right")
                        else:
                            try:
                                unc_value = float(child.text)
                            except TypeError: # empty tag is of type NULL
                                self.warning('<uncertainty> tag is empty; ' +
                                             'setting to 0')
                                unc_value = 0.
                            except ValueError:
                                raise ExpInputError(self.filepath,
                                                    "value of <uncertainty> tag is " +
                                                    "not a number")
                                
                        param["uncertainty"][child.attrib["side"]] = unc_value
                else:
                    raise ExpInputError(self.filepath,
                                        "subtag or param should be uncertainty")
                    
            elif dim == 2 and type == "n":
                allowed_tags = ["a", "b", "c"]
                if child.tag not in allowed_tags:
                    raise ExpInputError(self.filepath,
                                        "only allowed tags are <a>, <b> and <c> in " +
                                        "block param in 2D normal mode")

                if child.tag in param:
                    self.warning("redefinition of tag <" + child.tag + ">")
                    
                try:
                    param_value = float(child.text)
                except TypeError: # empty tag is of type NULL
                    self.warning('<' + child.tag + '> tag is empty; ' +
                                 'setting to 0')
                    param_value = 0.
                except ValueError:
                    raise ExpInputError(self.filepath,
                                        "value of <" + child.tag + "> tag is not a number")
                
                param[child.tag] = param_value

            elif dim == 2 and (type == "vn" or type == "p"):
                allowed_tags = ["uncertainty", "correlation"]
                if child.tag not in allowed_tags:
                    raise ExpInputError(self.filepath,
                                        "only allowed tags are <uncertainty> and <correlation> in " +
                                        "block param in 2D variable normal mode")

                if child.tag == "uncertainty":
                    if child.attrib["axis"] not in ["x", "y"]:
                        raise ExpInputError(self.filepath,
                                         "axis attribute of uncertainty is not x nor y")
                    elif child.attrib["side"] not in ["left", "right"]:
                        raise ExpInputError(self.filepath,
                                         "side attribute of uncertainty is not left nor right")
                    else:
                        try:
                            unc_value = float(child.text)
                        except TypeError: # empty tag is of type NULL
                            self.warning('<uncertainty> tag is empty; ' +
                                         'setting to 0')
                            unc_value = 0.
                        except ValueError:
                            raise ExpInputError(self.filepath,
                                         "value of <uncertainty> tag is not a number")

                        axis_label = child.attrib["axis"]
                        side_label = child.attrib["side"]
                        param[child.tag][axis_label][side_label] = unc_value

                if child.tag == "correlation":
                    try:
                        corr_value = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        self.warning('<correlation> tag is empty; ' +
                                         'setting to 0')
                        corr_value = 0.
                    except ValueError:
                        raise ExpInputError(self.filepath,
                                         "value of <correlation> tag is not a number")

                    param[child.tag] = corr_value
            elif dim >= 3:
                allowed_tags = ["uncertainty", "correlation"]
                if child.tag not in allowed_tags:
                    raise ExpInputError(self.filepath,
                                        "only allowed tags are <uncertainty> and <correlation> in " +
                                        "block param in 2D variable normal mode")

                if child.tag == "uncertainty":
                    if child.attrib["axis"] not in axis_allowedtags:
                        raise ExpInputError(self.filepath,
                                         "axis attribute of uncertainty is not correct")
                    elif child.attrib["side"] not in ["left", "right"]:
                        raise ExpInputError(self.filepath,
                                         "side attribute of uncertainty is not left nor right")
                    else:
                        try:
                            unc_value = float(child.text)
                        except TypeError: # empty tag is of type NULL
                            self.warning('<uncertainty> tag is empty; ' +
                                         'setting to 0')
                            unc_value = 0.
                        except ValueError:
                            raise ExpInputError(self.filepath,
                                         "value of <uncertainty> tag is not a number")

                        axis_label = child.attrib["axis"]
                        side_label = child.attrib["side"]
                        param[child.tag][axis_label][side_label] = unc_value

                if child.tag == "correlation":
                    if child.attrib["entry"] not in corr_allowedtags:
                        raise ExpInputError(self.filepath,
                                         "entry attribute of correlation is not correct")
                    else:
                        try:
                            corr_value = float(child.text)
                        except TypeError: # empty tag is of type NULL
                            self.warning('<correlation> tag is empty; ' +
                                         'setting to 0')
                            corr_value = 0.
                        except ValueError:
                            raise ExpInputError(self.filepath,
                                         "value of <correlation> tag is not a number")

                        entry_label = child.attrib["entry"]
                        param[child.tag][entry_label] = corr_value

        if dim == 1 and type == "p":
            sigm = abs(param["uncertainty"]["left"])
            sigp = param["uncertainty"]["right"]
            param["gamma"] = solve_bifurcation_f_gamma(sigm,sigp,1000)
            param["nu"] = 0.5/(param["gamma"]*sigp - np.log(1+param["gamma"]*sigp))

        if dim == 2 and type == "p":
            p = param["correlation"]
            sig1p = param["uncertainty"]["x"]["right"]
            sig1m = abs(param["uncertainty"]["x"]["left"])
            sig2p = param["uncertainty"]["y"]["right"]
            sig2m = abs(param["uncertainty"]["y"]["left"])

            param["gamma"]["x"] = solve_bifurcation_f_gamma(sig1m,sig1p,1000)
            param["gamma"]["y"] = solve_bifurcation_f_gamma(sig2m,sig2p,1000)

            param["nu"]["x"] = 0.5/(param["gamma"]["x"]*sig1p - np.log(1+param["gamma"]["x"]*sig1p))
            param["nu"]["y"] = 0.5/(param["gamma"]["y"]*sig2p - np.log(1+param["gamma"]["y"]*sig2p))

            param["A_corr"] = fsolve(f_Poisson_corr,0,args=(p,param["nu"]["x"],param["nu"]["y"]))[0]
            param["alpha_corr"] = np.log(1+param["A_corr"])

        if dim >= 3:
        # define the correlation matrix:
            corr_m = np.array([[1., param["correlation"]["d1d2"], param["correlation"]["d1d3"]],
                               [param["correlation"]["d1d2"],  1, param["correlation"]["d2d3"]],
                               [param["correlation"]["d1d3"], param["correlation"]["d2d3"],  1]])
            unc_right = np.array([[param["uncertainty"]["d1"]["right"],param["uncertainty"]["d2"]["right"],param["uncertainty"]["d3"]["right"]]])
            unc_left = np.array([[param["uncertainty"]["d1"]["left"],param["uncertainty"]["d2"]["left"],param["uncertainty"]["d3"]["left"]]])
            for i in range(4,dim+1):
                unc_right = np.append(unc_right,[[param["uncertainty"][axis_allowedtags[i-1]]["right"]]],axis=1)
                unc_left = np.append(unc_left,[[param["uncertainty"][axis_allowedtags[i-1]]["left"]]],axis=1)
                col = np.array([[]])
                for j in range(1,i):
                    r = "d" + str(j)
                    col = np.append(col,[[param["correlation"][r+axis_allowedtags[i-1]]]],axis=1)
                row = np.append(col,[[1.]],axis=1) 
                corr_m = np.concatenate((corr_m,col.T),axis=1)
                corr_m = np.concatenate((corr_m,row),axis=0)

            if type == "n":
                unc_sym = (unc_right + abs(unc_left))/2  # need symmetric errors
                cov_m = unc_sym*corr_m*unc_sym.T
                inv_cov_m = np.linalg.inv(cov_m)
                param["inv_cov_m"] = inv_cov_m
            elif type == "vn":
                param["VGau"] = unc_right*abs(unc_left)
                param["VGau_prime"] = unc_right - abs(unc_left)
                param["corr_m"] = corr_m

        # check that everything is there
        if (type == "n" or type == "vn" or type == "p") and dim == 1:
            if ("uncertainty" not in param or
                "left" not in param["uncertainty"] or
                "right" not in param["uncertainty"]):
                raise ExpInputError(self.filepath,
                                    "uncertainties are not given consistently in block param")
            if (param["uncertainty"]["left"] == 0 and param["uncertainty"]["right"] == 0):
                raise ExpInputError(self.filepath,
                                    "uncertainties are all zero")

        elif type == "n" and dim == 2:
            if "a" not in param or "b" not in param or "c" not in param:
                raise ExpInputError(self.filepath,
                                    "a, b, c tags are not given in block param")

        elif (type == "vn" or type == "p") and dim == 2:  # added by LDN
            if ("uncertainty" not in param or
                "x" not in param["uncertainty"] or
                "y" not in param["uncertainty"] or 
                "left" not in param["uncertainty"]["x"] or 
                "right" not in param["uncertainty"]["x"] or
                "left" not in param["uncertainty"]["y"] or 
                "right" not in param["uncertainty"]["y"] or
                "correlation" not in param):
                raise ExpInputError(self.filepath,
                             "uncertainty or correlation tags are not given correctly in block param")
            if (param["uncertainty"]["x"]["left"] == 0 and param["uncertainty"]["x"]["right"] == 0):
                raise ExpInputError(self.filepath,
                                    "x uncertainties are all zero")
            if (param["uncertainty"]["y"]["left"] == 0 and param["uncertainty"]["y"]["right"] == 0):
                raise ExpInputError(self.filepath,
                                    "y uncertainties are all zero")
            if (abs(param["correlation"]) == 1):
                raise ExpInputError(self.filepath,
                                    "correlation is (minus) unity, cannot handle this case")
        elif dim >= 3:  # added by LDN
            if ("uncertainty" not in param or
                "correlation" not in param):
                raise ExpInputError(self.filepath,
                             "uncertainty or correlation tags are not given in block param")

        # or the grid
        grid = {}
        Lxy = None
        
        for child in root:
            if child == "grid":
                break
        grid_raw = child.text

        if type == "f" and dim == 1:
            x = []
            L = []
        
            grid_raw = grid_raw.strip("\n").strip().split("\n")
        
            i = -1
            for line in grid_raw:
                tab = line.split()
                if len(tab) != 2:
                    raise ExpInputError(self.filepath,
                                        'incorrect <grid> entry on line "' + line + '"')

                cur_x = float(tab[0])
                cur_L = float(tab[1])

                if cur_x not in x:
                    x.append(cur_x)
                    L.append(cur_L)
                    i += 1
                else:
                    i = x.index(cur_x)
            
            grid["x"] = x
            grid["L"] = L
            LChi2min = min(grid["L"])

            Lxy = interpolate.UnivariateSpline(grid["x"], grid["L"], k=3, s=0)

        elif type == "f" and dim == 2:
            x = []
            y = []
            L = []
            
            grid_raw = grid_raw.strip("\n").strip().split("\n")
            
            i = -1
            
            for line in grid_raw:
                tab = line.split()
                if len(tab) != 3:
                    raise ExpInputError(self.filepath,
                                        'incorrect <grid> entry on line "' + line + '"')

                cur_x = float(tab[0])
                cur_y = float(tab[1])
                cur_L = float(tab[2])

                if cur_x not in x:
                    x.append(cur_x)
                    L.append([])
                    i += 1
                else:
                    i = x.index(cur_x)

                if cur_y not in y:
                    y.append(cur_y)

                L[i].append(cur_L)

            grid["x"] = np.array(x)
            grid["y"] = np.array(y)
            grid["L"] = np.array(L)
            
            LChi2min = min(min(p[1:]) for p in grid["L"])

            Lxy = interpolate.RectBivariateSpline(grid["x"],
                    grid["y"], grid["L"])
            
        return (bestfit, param, grid, Lxy, LChi2min)

    def check_multiprod(self, eff_dict, multiprod):
        """..."""

        # check consistency in the definition of multi-particle labels
        for (prod,decay) in eff_dict:
            if prod in multiprod:
                # there is a multi-particle label
                # in that case, check if individual particle are also defined
                for label in multiprod[prod]:
                    if (label,decay) in eff_dict:
                        raise ExpInputError(self.filepath,
                                            '<eff> tags for "' + label + '" and "' +
                                            prod + '" cannot both be defined')
                # also, overlapping multi-particle labels cannot be both used (VH or VVH),
                # not both
                if prod == "ZH" and ("VH",decay) in eff_dict:
                    raise ExpInputError(self.filepath,
                                            '<eff> tags for "' + "VH" + '" and "' +
                                             prod + '" cannot both be defined')
                if prod == "ZH" and ("VVH",decay) in eff_dict:
                    raise ExpInputError(self.filepath,
                                            '<eff> tags for "' + "VVH" + '" and "' +
                                             prod + '" cannot both be defined')
                if prod == "VH" and ("VVH",decay) in eff_dict:
                    raise ExpInputError(self.filepath,
                                            '<eff> tags for "' + "VVH" + '" and "' +
                                             prod + '" cannot both be defined')
                if prod == "tH" and ("top",decay) in eff_dict:
                    raise ExpInputError(self.filepath,
                                            '<eff> tags for "' + "top" + '" and "' +
                                             prod + '" cannot both be defined')

        # it is consistent, resolve multi-particle labels
        new_eff = {}
        for (prod,decay) in eff_dict:
            if prod in multiprod:
                for label in multiprod[prod]:
                    new_eff[label,decay] = eff_dict[prod,decay]*multiprod[prod][label]

        for elem in new_eff:
            eff_dict[elem] = new_eff[elem]


#-- functions for variable Gaussian and Poission likelihoods --

def solve_bifurcation_f_gamma(m, p, N):
    a = 0.
    b = 1.0/m
    for i in range(N):
        x = (a+b)/2 
        if np.exp(-x*(m+p)) <= (1-m*x)/(1+p*x): 
            a = x
        else:
            b = x
    return a

def f_Poisson_corr(y, *params):
    corr, z1, z2 = params
    f = z1*z2*y - corr*np.sqrt(z1*z2*(1+z2*(np.exp(z1*y**2) - 1)))
    return f

