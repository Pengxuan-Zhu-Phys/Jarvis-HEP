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
import os
try:
    from lxml import etree
except:
    import xml.etree.ElementTree as etree
from ..errors import UserInputError, HiggsMassError
from warnings import warn

class ReadUserInput:
    """Read the XML input and extracts all information."""

    modes = ["signalstrengths", "reducedcouplings"]

    def __init__(self, inputstring):
        """Initialize the reading of the user input from the XML input contained
        in the string inputstring."""

        self.root = etree.fromstring(inputstring)
    
        if self.root.tag != "lilithinput":
            raise UserInputError('root tag is not <lilithinput>')

        self.mode, n_higgses = self.getmode()

        # self.redC is for reduced coupling mode only
        # it is a list of Higgs particles containing the reduced
        # couplings and the extraBR as well as precision
        self.redC = []
        # self.mu is for signal strengths mode only
        # it is a list of Higgs particles containing the signal strengths
        # for the various (prod,decay) combinations 
        self.mu = []

        for i in range(n_higgses):
            if self.mode == "reducedcouplings":
                self.redC.append(self.get_nextreducedcouplings())
            elif self.mode == "signalstrengths":
                self.mu.append(self.get_nextsignalstrengths())

    def warning(self, message):
        """Customized warnings."""

        warn(message, Warning, stacklevel=3)

    def getmode(self):
        """Get the mode (signal strengths or reduced couplings) from the
           XML tree."""

        curmode = ""

        num_higgses = 0
        name_higgses = []

        for child in self.root:
            if child.tag is etree.Comment:
                # ignore all comments
                continue
            if child.tag not in ReadUserInput.modes:
                # if there is some other tag
                self.warning('<' + child.tag + '> is present outside of a <' +
                             ReadUserInput.modes[0] + '> or <' +
                             ReadUserInput.modes[1] + '> block and will be ' +
                             'ignored')
            else:
                num_higgses += 1
                if curmode == "":
                    curmode = child.tag
                    # next, handle possible multiparticle input
                    if "part" in child.attrib and child.attrib["part"] != "":
                       name_higgses.append(child.attrib["part"])
                else:
                    # if more than one mode tag is defined in the input
                    if child.tag != curmode:
                        raise UserInputError(
                            '<' + child.tag + '> and <' + curmode +
                            '> blocks cannot both be present')
                    else:
                        if ("part" in child.attrib and
                            child.attrib["part"] != ""):
                            name_higgses.append(child.attrib["part"])
    
        # check if at least one mode is present
        if curmode == "":
            raise UserInputError(
                'no mode block has been defined (' +
                '<signalstrengths> or <reducedcouplings>)')

        # check if two higgses have the same name
        if any(name_higgses.count(x) > 1 for x in name_higgses):
            self.warning('there are several higgses with the same name; ' +
                         'we will consider that all particles are different ' +
                         'and contribute to the signal')

        return curmode, num_higgses
    
    def getmass(self, higgs_block):
        """Get the Higgs mass (in GeV) from the XML tree."""
        
        mass = -1.
        default_mass = 125.09
        min_mass = 123.
        max_mass = 128.

        for child in higgs_block: # <reducedcouplings> or <signalstrengths> tag
            if child.tag == "mass":
                if mass > 0:
                    self.warning('redefinition of the Higgs mass for the ' +
                                 'particle ' + higgs_block.tag)
                try:
                    mass = float(child.text)
                    assert mass >= min_mass and mass <= max_mass
                except TypeError: # empty tag is of type NULL
                    mass = default_mass
                    self.warning('<mass> tag is empty; setting to ' +
                                 'SM value=' + default_mass + ' GeV')
                except ValueError:
                    raise UserInputError(
                        'value of the <mass> tag is not a number.')
                except AssertionError:
                    raise HiggsMassError('<mass> is not between ' + str(min_mass) +
                                         ' and ' + str(max_mass) + ' GeV.')

        # if no <mass> block, set to default value
        if mass < 0:
            mass = default_mass

        return mass

    def get_nextreducedcouplings(self):
        """..."""

        redCp = {"extra": {"precision": ""}}
        
        accepted_C = ["ff", "uu", "ll", "tt", "cc", "dd", "bb", "tautau", "VV", "WH", "ZH",
                      "WW", "ZZ", "gg", "gammagamma", "Zgamma", "mumu", "VBF"]
        accepted_precision = ["LO", "BEST-QCD"]
        accepted_extraBR = ["invisible", "undetected"]

        red_coupl = None

        n = 0
        for child in self.root:
            if child.tag != "reducedcouplings":
                continue
            if n == len(self.redC):
                red_coupl = child
                if "part" in child.attrib and child.attrib["part"] != "":
                    redCp["extra"]["name"] = child.attrib["part"]
            n += 1

        if red_coupl is None:
            # this is not supposed to happen
            raise UserInputError(
                'no <reducedcouplings> tag is matching the requested ' +
                 'Higgs particle')

        redCp["extra"]["mass"] = self.getmass(red_coupl)

        for child in red_coupl:
            if child.tag == "C":
                if "to" not in child.attrib:
                    self.warning('attribute "to" is missing in <C> tag')
                elif child.attrib["to"] not in accepted_C:
                    self.warning('<C> tag to "' + child.attrib["to"] +
                                 '" is unknown')
                else:
                    if child.attrib["to"] == "gg":
                        # gg has optional attribute for="prod|decay|all"
                        if ("for" not in child.attrib or
                            child.attrib["for"] == "all"):
                            C_to = "gg_all"
                        elif child.attrib["for"] == "decay":
                            C_to = "gg_decay"
                        elif child.attrib["for"] == "prod":
                            # gg for="prod" has optional attribute at=lhc8|all
                            if ("at" not in child.attrib or
                                child.attrib["at"] == "all"):
                                C_to = "gg_prod_all"
                            elif child.attrib["at"] == "lhc8":
                                C_to = "gg_prod_lhc8"
                            else:
                                raise UserInputError(
                                    '<C to="gg" for="prod">' +
                                    ' and at="' + child.attrib["at"] +
                                    '" is unknown...')
                        else:
                            raise UserInputError(
                                '<C to="gg"> for "' + child.attrib["for"] +
                                '" is unknown...')
                    else:
                        C_to = child.attrib["to"]
                    
                    if child.attrib["to"] in ["ff", "uu", "ll", "tt", "cc", "dd", "bb", "tautau","mumu"]:
                        # if fermion, there could be a real and imaginary part
                        if "part" in child.attrib:
                            if child.attrib["part"] in ["re", "im"]:
                                C_to += "_" + child.attrib["part"]
                    
                    if C_to in redCp:
                        self.warning('reduced coupling for "' +
                        C_to + '" is being redefined')
                    try:
                        if C_to == "gg_all":
                            redCp["gg_decay"] = float(child.text)
                            redCp["gg_prod_lhc8"] = float(child.text)
                            redCp["gg_prod_lhc13"] = float(child.text)
                        elif C_to == "gg_prod_all":
                            redCp["gg_prod_lhc8"] = float(child.text)
                        else:
                            redCp[C_to] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        redCp[C_to] = 1.
                        self.warning('<C> tag to "' + child.attrib["to"] +
                                '" is empty; setting to SM value = 1')
                    except ValueError:
                        raise UserInputError(
                            'value of the <C> tag to "' + child.attrib["to"] +
                            '" is not a number.')
        
            if child.tag == "precision":
                if child.text not in accepted_precision:
                    self.warning('"' + str(child.text) +
                            '" precision is not allowed; setting to default ' +
                            'value = "BEST-QCD"')
                    redCp["extra"]["precision"] = "BEST-QCD"
                else:
                    if redCp["extra"]["precision"] != "":
                        self.warning('precision is being redefined')
                    redCp["extra"]["precision"] = child.text

            if child.tag == "extraBR":
                for subchild in child:
                    if "to" not in subchild.attrib:
                        self.warning('attribute "to" is missing in <BR> tag')
                    elif subchild.attrib["to"] not in accepted_extraBR:
                        self.warning('<BR> tag for type "' + subchild.attrib["to"] +
                                     '" is unknown.')
                    else:
                        if "BR" + subchild.attrib["to"] in redCp["extra"]:
                            self.warning('<BR> tag with type ' + subchild.attrib["to"] +
                                         ' is being redefined')
                        
                        try:
                            redCp["extra"]["BR" + subchild.attrib["to"]] = float(subchild.text)
                        except TypeError: # empty tag is of type NULL
                            self.warning('<BR> tag to "' + subchild.attrib["to"] +
                                         '" is empty; setting to SM value = 0')
                            redCp["extra"]["BR" + subchild.attrib["to"]] = 0.
                        except ValueError:
                            raise UserInputError(
                                'value of the <BR> tag for "' +
                                subchild.attrib["to"] + '" is not a ' +
                                'number.')

        # --- putting together real and imaginary part when present
        new_redC = {}
        original_redC = {}
        for p in redCp:
            new_p = 0.
            if p[-3:] in ["_im", "_re"]:
                if p[:-3] in redCp:
                    raise UserInputError(
                        p[:-3] + ' and ' + p + ' are defined at the same time')
                
                if p[-3:] == "_re":
                    num = redCp[p]
                else:
                    num = redCp[p]*1j
                
                if p[:-3] in new_redC:
                    new_redC[p[:-3]] += num
                else:
                    new_redC[p[:-3]] = num
            else:
                original_redC[p] = redCp[p]

        redCp = original_redC
        redCp.update(new_redC)

        # ----------------------------------
        # checking consistency of the input
        # ----------------------------------

        # --- set extraBR invisible and undetected to 0 if not given
        if "BRinvisible" not in redCp["extra"]:
            redCp["extra"]["BRinvisible"] = 0.
        if "BRundetected" not in redCp["extra"]:
            redCp["extra"]["BRundetected"] = 0.

        # --- checking the multiparticles
        multiparticles = {"VV": ["WW", "ZZ"], "uu": ["tt", "cc"],
                          "dd": ["bb", "tautau","mumu"], "ll": ["tautau","mumu"]}
        multiparticles2 = {"ff": ["tt", "cc", "bb", "tautau","mumu"]}
        # two separate dicts because "ff" is a special case since it overlaps
        # with the multiparticle labels "uu", "dd" and "ll"---it should be treated
        # in the very end

        for multip, p_list in list(multiparticles.items()):
            self.check_multiparticle(redCp, multip, p_list)
        
        for multip, p_list in list(multiparticles2.items()):
            self.check_multiparticle(redCp, multip, p_list)

        # now all reduced couplings have been properly defined, one can
        # delete all multiparticle labels
        redCclean = redCp.copy()
        for p in redCp:
            if p in multiparticles or p in multiparticles2:
                del redCclean[p]
        redCp = redCclean

        # --- checking that the couplings to all particles are well defined in
        #     input (not including optionnal couplings for the the loop-induced
        #     processes)
        mandatory_particles = ["WW", "ZZ", "tt", "cc", "mumu", "bb", "tautau", "qqZH", "WH"]

        for p in mandatory_particles:
            if p not in redCp:
                #self.warning('reduced coupling to "' + p + '" is not specified; fixing it to SM value=1')
                if p=="qqZH":
                    redCp[p]=redCp["ZZ"]
                elif p == "WH":
                    redCp[p] = redCp["WW"]
                else:
                    redCp[p] = 1.

        # --- checking precision
        if redCp["extra"]["precision"] == "":
            #self.warning('precision has not been specified; setting to default value = "BEST-QCD"')
            redCp["extra"]["precision"] = "BEST-QCD"

        # --- setting C_VBF = C_W in case C_W = C_Z
        if "VBF" not in redCp:
            if abs(redCp["WW"] - redCp["ZZ"]) < 1e-6:
                redCp["VBF"] = redCp["WW"]

        return redCp

    def get_nextsignalstrengths(self):
        """..."""

        mup = {"extra": {}}

        accepted_prod = ["ggH", "VVH", "ttH", "VBF", "VH", "WH", "qqZH", "ggZH", "ZH", "tHq", "tHW", "tH", "top", "bbH"]
        accepted_decay = ["gammagamma", "VV", "WW", "ZZ", "bb", "tautau",
                          "dd", "uu", "ll", "cc", "ff", "Zgamma", "mumu", "invisible", "gg"]

        mu_block = None

        n = 0
        for child in self.root:
            if child.tag != "signalstrengths":
                continue
            if n == len(self.mu):
                mu_block = child
                if "part" in child.attrib and child.attrib["part"] != "":
                    mup["extra"]["name"] = child.attrib["part"]
            n += 1

        mup["extra"]["mass"] = self.getmass(mu_block)

        for child in mu_block:
            if child.tag == "mu":
                if "prod" not in child.attrib or "decay" not in child.attrib:
                    self.warning('attribute "prod" or "decay" is ' +
                                 'missing in <mu> tag')
                elif child.attrib["prod"] not in accepted_prod:
                    self.warning('<mu> tag with prod="' +
                                 child.attrib["prod"] + '" and decay="' +
                                 child.attrib["decay"] + '" has unknown prod')
                elif child.attrib["decay"] not in accepted_decay:
                    self.warning('<mu> tag with prod="' +
                                 child.attrib["prod"] + '" and decay="' +
                                 child.attrib["decay"] + '" has unknown decay')
                else:
                    prod = child.attrib["prod"]
                    decay = child.attrib["decay"]
                    if (prod,decay) in mup:
                        self.warning('<mu> tag with prod="' + prod +
                                     '" and decay="' + decay + '" is being ' +
                                     'redefined')

                    try:
                        mup[prod,decay] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        mup[prod, decay] = 1.
                        self.warning('<mu> tag with prod="' +
                                     prod + '" and decay="' +
                                     decay + '" is empty; ' +
                                     'setting to SM value=1')
                    except ValueError:
                        raise UserInputError(
                            'value of the <mu> tag with prod="' + prod +
                            '" and decay="' + decay + '" is not a number.')
            elif child.tag == "redxsBR":
                if "prod" not in child.attrib or "decay" not in child.attrib:
                    self.warning('attribute "prod" or "decay" is ' +
                                 'missing in <redxsBR> tag')
                elif child.attrib["prod"] not in accepted_prod:
                    self.warning('<mu> tag with prod="' +
                                 child.attrib["prod"] + '" and decay="' +
                                 child.attrib["decay"] + '" has unknown prod')
                elif child.attrib["decay"] != "invisible":
                    self.warning('<redxsBR> tag with prod="' +
                                 child.attrib["prod"] + '" and decay="' +
                                 child.attrib["decay"] + '" has unknown decay')
                else:
                    prod = child.attrib["prod"]
                    decay = child.attrib["decay"]

                    if (prod,decay) in mup:
                        self.warning('<redxsBR> tag with prod="' + prod +
                                     '" and decay="' + decay + '" is being ' +
                                     'redefined')

                    try:
                        mup[prod,decay] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        mup[prod, decay] = 1.
                        self.warning('<redxsBR> tag with prod="' +
                                 prod + '" and decay="' +
                                 decay + '" is empty; ' +
                                 'setting to SM value=1')
                    except ValueError:
                        raise UserInputError(
                            'value of the <mu> tag with prod="' + prod +
                            '" and decay="' + decay + '" is not a number.')
        
        # ----------------------------------
        # checking consistency of the input
        # ----------------------------------

        # --- checking the multiparticle labels
        multiprod = {"VH": ["WH", "qqZH", "ggZH"]}
        multiprod2 = {"VVH": ["VBF", "WH", "qqZH", "ggZH"]}
        multiprod3 = {"tH": ["tHq", "tHW"]}
        multiprod4 = {"top": ["ttH", "tHq", "tHW"]}
        multiprod5 = {"ZH": ["qqZH", "ggZH"]}
        multidecay = {"VV": ["WW", "ZZ"], "uu": ["cc"], "dd": ["bb", "tautau","mumu"], "ll": ["tautau","mumu"]}
        multidecay2 = {"ff": ["cc", "bb", "tautau","mumu"]}
        # two separate dicts for prod and decay because "ff" is a special case
        # since it overlaps with the multiparticle labels "uu", "dd" and "ll"
        # also "VVH" includes "VH"; this should be treated at the end

        # -- first: production
        for multip, p_list in list(multiprod.items()):
            self.check_multiprod(mup, multip, p_list)

        for multip, p_list in list(multiprod2.items()):
            self.check_multiprod(mup, multip, p_list)

        for multip, p_list in list(multiprod3.items()):
            self.check_multiprod(mup, multip, p_list)

        for multip, p_list in list(multiprod4.items()):
            self.check_multiprod(mup, multip, p_list)

        for multip, p_list in list(multiprod5.items()):
            self.check_multiprod(mup, multip, p_list)

        # now, one can delete all multiprod labels
        for key,mu_value in list(mup.items()):
            if key == "extra":
                continue
            prod,decay = key
            if prod in multiprod or prod in multiprod2 or prod in multiprod3 or prod in multiprod4 or prod in multiprod5 :
                del mup[prod,decay]

        # -- then: decay
        for multip, p_list in list(multidecay.items()):
            self.check_multidecay(mup, multip, p_list)

        for multip, p_list in list(multidecay2.items()):
            self.check_multidecay(mup, multip, p_list)

        # now, one can delete all multidecay labels
        for key,mu_value in list(mup.items()):
            if key == "extra":
                continue
            prod,decay = key
            if decay in multidecay or decay in multidecay2:
                del mup[prod,decay]

        # --- checking that all mandatory signal strengths are present in the
        #     input
        mandatory_prod = ["ggH", "VBF", "WH", "qqZH", "ggZH", "ttH", "tHq", "tHW", "bbH"]
        mandatory_decay = ["gammagamma", "ZZ", "WW", "bb", "tautau", "gg", "cc", "mumu", "Zgamma", "invisible"]

        warning_prod = ["ggH", "VBF", "WH", "qqZH", "ttH"]
        warning_decay = ["gammagamma", "ZZ", "WW", "bb", "tautau"]

        mandatory_mus = []
        for prod in mandatory_prod:
            for decay in mandatory_decay:
                mandatory_mus.append((prod,decay))

        for (prod,decay) in mandatory_mus:
            if (prod,decay) not in mup:
                if decay == "invisible":
                    mup[prod,decay] = 0.
                else:
                    mup[prod,decay] = 1.
                if prod in warning_prod and decay in warning_decay:
                    self.warning('signal strength for prod "' + prod +
                             '" and decay "' + decay + '" is not  specified; ' +
                             'fixing it to SM value=1')

        return mup

    def check_multiparticle(self, redCp, multip, p_list):
        """for reduced couplings"""
        
        if multip in redCp:
            if set(p_list).issubset(redCp):
                # if multiparticle and individual particles are all
                # defined, check consistency
                for p in p_list:
                    if redCp[multip] != redCp[p]:
                        self.warning('inconsistent definition of the ' +
                                     'couplings to "' + multip + '" and "' +
                                     p + '"; skipping <C to="' +
                                     multip + '"> tag')
            else:
                # only multiparticle label is present, or multiparticle
                # label and tags for only part of the associated particles
                given_p = [part for part in p_list if part in redCp]
                
                # check if a particle tag is given inconsistently
                for p in given_p:
                    if redCp[multip] != redCp[p]:
                        raise UserInputError(
                            'inconsistent definition of the ' +
                            'couplings to "' + multip + '" and to "' +
                             p + '"')
                
                if len(given_p) == 0:
                    # only multiparticle tag is defined
                    for p in p_list:
                        redCp[p] = redCp[multip]
                else:
                    # when multiparticle tag and part of the particle tags
                    # are given, with equal values
                    other_p = [part for part in p_list if
                               part not in redCp]
                    for op in other_p:
                        redCp[op] = redCp[multip]
                    self.warning('couplings to ' + multip + ' and ' +
                                 str(given_p) + ' are both defined and ' +
                                 'equal; assuming that the couplings to ' +
                                 str(other_p) + ' are the same')

    def check_multiprod(self, mup, multip, p_list):
        """for signal strengths mode"""

        for key,mu_value in list(mup.items()):
            if key == "extra":
                continue
            prod,decay = key
            if prod == multip:
                if set([(subprod,decay) for subprod in p_list]).issubset(mup):
                    # if multiprod and individual prod are all
                    # defined, check consistency
                    for subprod in p_list:
                        if mu_value != mup[subprod,decay]:
                            self.warning('inconsistent definition of the ' +
                                         'mu prod="' + multip + '" and "' +
                                         prod + '" when decay="' + decay + '"; ' +
                                         'skipping ' + 'mu prod="' + multip + '" tag')
                else:
                    # only multiprod label is present, or multiprod
                    # label and tags for only part of the associated prod
                    given_mu = [(subprod, decay) for subprod in p_list if
                                (subprod, decay) in mup]
                    
                    # check if a particle tag is given inconsistently
                    for (subprod,decay) in given_mu:
                        if mu_value != mup[subprod,decay]:
                            raise UserInputError(
                                'inconsistent definition of the mu for prod "' +
                                subprod + '" and "' + multip +
                                '" when decay is "' + decay + '"')
            
                    if len(given_mu) == 0:
                        # only multiparticle tag is defined
                        for subprod in p_list:
                            mup[subprod,decay] = mup[multip,decay]
                    else:
                        # when multiprod tag and part of the prod tags
                        # are given, with equal values
                        other_mu = [(subprod, decay) for subprod in p_list if
                                    (subprod, decay) not in mup]
                        for (op,decay) in other_mu:
                            mup[op,decay] = mup[multip,decay]
                        self.warning('signal strengths with prod ' + multip +
                                     ' and ' +
                                     str(given_mu) + ' are all defined and ' +
                                     'equal; assuming that the signal strengths to '+
                                     str(other_mu) + ' are the same [for decay = "' +
                                     decay + '"]')

    def check_multidecay(self, mup, multip, p_list):
        """for signal strengths mode"""

        for key,mu_value in list(mup.items()):
            if key == "extra":
                continue
            prod,decay = key
            if decay == multip:
                if set([(prod,subdecay) for subdecay in p_list]).issubset(mup):
                    # if multidecay and individual decay are all
                    # defined, check consistency
                    for subdecay in p_list:
                        if mu_value != mup[prod,subdecay]:
                            self.warning('inconsistent definition of the ' +
                            'mu decay="' + multip + '" and "' +
                            decay + '" when prod="' + prod + '"; ' +
                            'skipping ' + 'mu decay="' + multip + '" tag')
                else:
                    # only multidecay label is present, or multidecay
                    # label and tags for only part of the associated prod
                    given_mu = [(prod, subdecay) for subdecay in p_list if
                                (prod, subdecay) in mup]
                    
                    # check if a particle tag is given inconsistently
                    for (prod,subdecay) in given_mu:
                        if mu_value != mup[prod,subdecay]:
                            raise UserInputError(
                                'inconsistent definition of mu for decay "' +
                                subdecay + '" and "' + multip +
                                '" when prod is "' + prod + '"')
            
                    if len(given_mu) == 0:
                        # only multiparticle tag is defined
                        for subdecay in p_list:
                            mup[prod,subdecay] = mup[prod,multip]
                    else:
                        # when multidecay tag and part of the decay tags
                        # are given, with equal values
                        other_mu = [(prod, subdecay) for subdecay in p_list if
                                    (prod, subdecay) not in mup]
                        for (prod,od) in other_mu:
                            mup[prod,od] = mup[prod,multip]
                        self.warning('signal strengths with decay ' + multip +
                                     ' and ' +
                                     str(given_mu) + ' are all defined and ' +
                                     'equal; assuming that the signal strengths to '+
                                     str(other_mu) + ' are the same [for prod = "' +
                                     prod + '"]')

