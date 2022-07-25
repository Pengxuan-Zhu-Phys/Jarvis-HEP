#! /usr/bin/env python

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

# standard modules
import os.path, time, sys
import importlib
from warnings import warn
# Lilith library
from .errors import ExpNdfComputationError, UserMuTotComputationError, \
                   UserInputIOError
from .internal.readexpinput import ReadExpInput
from .internal.readuserinput import ReadUserInput
from .internal.computereducedcouplings import ComputeReducedCouplings
from .internal.computemufromreducedcouplings import \
    ComputeMuFromReducedCouplings
from .internal.computelikelihood import compute_likelihood
import lilith.internal.writeoutput as writeoutput
import lilith.version as version

class Lilith:
    """Main class. Reads the experimental and user input and computes the
    likelihood.
    """

    default_exp_list = ("/".join(
        os.path.dirname(os.path.abspath(__file__)).split("/")[:-1]) +
        "/data/latest.list")

    def __init__(self, verbose=False, timer=False):
        """Initialize the relevant attributes."""

        # controls the information displayed on the screen
        self.verbose = verbose
        self.timer = timer

        # objects needed for the computation of reduced couplings and of
        # signal strengths from reduced couplings
        self.coupling_computation = None
        self.mu_computation = None

        # information read from the database of experimental results
        # each element of self.exp_mu corresponds to an XML file
        self.exp_mu = []
        self.exp_ndf = 0
        self.dbversion = "??.??"

        # information read from the user input
        # - in signal strengths mode, self.couplings remains empty
        # - each element of self.couplings and self.user_mu
        #   corresponds to a Higgs particle defined in the input
        self.mode = ""
        self.couplings = []
        self.user_mu = []
        self.user_mu_tot = {}

        # detailed likelihood results as well as the final result
        # self.l is defined as -2 log(total likelihood)
        self.results = []
        self.l = 0.
        self.l_SM = 0.

    def info(self, message):
        """Print information only is verbose is True"""
        
        if self.verbose:
            print(message)

    def tinfo(self, action, dt):
        """Print time taken for a given action"""

        if self.timer:
            print("- " + action + ": " + str(dt) + "s")

    def readuserinput(self, userinput):
        """Read the XML input given by the user."""

        self.info("Reading the user input...")
        t0 = time.time()
        userinput = ReadUserInput(userinput)
        self.tinfo("reading of the user input", time.time() - t0)
        self.mode = userinput.mode

        if userinput.mode == "reducedcouplings":
            self.user_mu = []
            self.user_mu_tot = {}
            self.couplings = userinput.redC
            self.info("User input: reduced couplings\n")
        else:
            self.couplings = []
            self.user_mu = userinput.mu
            self.compute_user_mu_tot()
            self.info("User input: signal strengths\n")

    def readuserinputfile(self, filepath):
        """Read the XML input given by the user in a file."""

        try:
            with open(filepath, "r") as f:
                self.readuserinput(f.read())
        except IOError as e:
            raise UserInputIOError(
                'I/O error({0}): {1}'.format(e.errno, e.strerror) + '; cannot' +
                ' open the user input file "' + filepath + '".')

    def computecouplings(self):
        """Computes missing reduced couplings."""

        computablecouplings = ["gg_prod_lhc8", "gg_decay", "gammagamma",
                               "Zgamma", "VBF","gg_prod_lhc13", "VBF13"]

        t0 = time.time()
        for n,redCp in enumerate(self.couplings, start=1):
            missing_couplings = False
            for coupling in computablecouplings:
                if coupling not in redCp:
                    missing_couplings = True
                    break

            if missing_couplings:
                # if reduced couplings need to be calculated for at least one of
                # the loop-induced processes or for VBF
                if self.coupling_computation is None:
                    self.coupling_computation = ComputeReducedCouplings(redCp)
                else:
                    self.coupling_computation.reset(redCp)

                new_redC = self.coupling_computation.getcouplings(redCp)

                info_str = ("The following reduced couplings have been "
                            " computed at " + redCp["extra"]["precision"] + 
                            " accuracy for the Higgs particle " + str(n))
                if "name" in redCp["extra"]:
                    info_str += " (" + redCp["extra"]["name"] + ")"
                info_str += ":"
                self.info(info_str)
                for cname,cvalue in list(new_redC.items()):
                    self.info(". " + cname + " = " + str(cvalue))
                redCp.update(new_redC)
        self.tinfo("computing missing reduced couplings", time.time() - t0)

    def computemufromreducedcouplings(self):
        """Computes signal strengths from reduced couplings."""

        t0 = time.time()
        self.user_mu = []
        for redCp in self.couplings:
            if self.mu_computation is None:
                self.mu_computation = ComputeMuFromReducedCouplings(
                    redCp["extra"]["mass"])
            else:
                self.mu_computation.reset(redCp["extra"]["mass"])
            self.user_mu.append(self.mu_computation.getmu(redCp))
        self.compute_user_mu_tot()
        self.tinfo("computing mu from reduced couplings",
                   time.time() - t0)

    def compute_user_mu_tot(self):
        """Adds up the signal strengths obtained from the user input."""

        if not self.user_mu:
            raise UserMuTotComputationError(
                "user_mu is empty, read signal strengths user input " +
                "or compute signal strengths from couplings first")

        self.user_mu_tot = {}
        for mup in self.user_mu:
            for key in mup:
                if key == "extra": # only consider signal strengths
                    continue
                if key in self.user_mu_tot:
                    self.user_mu_tot[key] += mup[key]
                else:
                    self.user_mu_tot[key] = mup[key]

    def readexpinput(self, filepath=default_exp_list):
        """Read the experimental input specified in a list file."""
    
        self.info("Processing the experimental input...")
        self.readdbversion()
        # initialize the reading of the experimental input
        exp_input = ReadExpInput()
        # read the list of XML files
        filelist = exp_input.get_filelist(filepath)
        # read and check each individual XML file
        for expfile in filelist:
            exp_input.read_file(expfile)
        self.exp_mu = exp_input.mu
        self.compute_exp_ndf()

    def readdbversion(self):
        dbversionfile = "/".join(os.path.dirname(os.path.abspath(__file__))
                          .split("/")[:-1])+"/data/version"
        try:
            with open(dbversionfile, "r") as f:
                self.dbversion = f.read().split("\n")[1]
        except IOError:
            warn("database version file data/version cannot be read", Warning)

    def compute_exp_ndf(self):
        self.exp_ndf = 0
        for mu in self.exp_mu:
            try:
                self.exp_ndf += mu["dim"]
            except KeyError:
                raise ExpNdfComputationError(
                    'there are missing elements in exp_mu: key ' + '"' + s +
                    '" is not found')

    def computelikelihood(self, userinput=None, exp_filepath=None,
                          userfilepath=None):
        """Computes the likelihood from the signal strengths (computed from)
           the user input and the experimental results."""

        if userinput is not None or userfilepath is not None:
            # read the user input and get user_mu_tot
            if userinput is not None:
                self.readuserinput(userinput)
            else:
                self.readuserinputfile(userfilepath)
            if self.couplings: # reduced coupling mode
                self.computecouplings()
                self.computemufromreducedcouplings()
        elif self.couplings and not self.user_mu_tot:
                self.computecouplings()
                self.computemufromreducedcouplings()

        if exp_filepath is not None:
            # read the experimental input and get exp_mu
            self.readexpinput(exp_filepath)
        elif not self.exp_mu:
            self.readexpinput()

        t0 = time.time()
        self.results, self.l = compute_likelihood(self.exp_mu,
                                                  self.user_mu_tot,self.mode)
        self.tinfo("computing the likelihood", time.time() - t0)
        
    def computeSMlikelihood(self, userinput=None, exp_filepath=None,
                          userfilepath=None):
        """Computes the SM likelihood from the signal strengths (computed from)
           the SM input and the experimental results."""
        
        self.readexpinput()
        decay_modes = ["gammagamma", "ZZ", "WW", "bb", "cc", "tautau", "Zgamma", "mumu", "gg","invisible"]
        prod_modes = ["ggH", "VBF", "WH", "qqZH", "ggZH", "ttH", "tHq", "tHW", "bbH"]
        SM_mu = dict(((l1,l2), float(l2!="invisible")) for l1 in prod_modes for l2 in decay_modes)
        self.results, self.l_SM = compute_likelihood(self.exp_mu, SM_mu, "signalstrengths")


    def writecouplings(self, filepath):
        writeoutput.couplings(self.couplings, filepath)

    def writesignalstrengths(self, filepath, tot=False):
        if tot:
            writeoutput.signalstrengths(self.user_mu_tot, filepath)
        else:
            writeoutput.signalstrengths(self.user_mu, filepath)

    def writeresults(self, filepath, slha=False):
        if slha:
            try:
                ndf = int(sys.argv[-2])
                l_ref = float(sys.argv[-1])
                if l_ref == -1:
                    l_ref = self.l_SM
                if l_ref == 0:
                    ndf = self.exp_ndf - ndf
                writeoutput.results_slha_pvalue(self.results, self.l, l_ref, ndf, filepath, self.dbversion)
# added by Ninh for v2.0
            except AttributeError:
                writeoutput.results_slha(self.results, self.l, self.l_SM, filepath)
            except IndexError:
                writeoutput.results_slha(self.results, self.l, self.l_SM, filepath)
# end of addition
            except ValueError:
                writeoutput.results_slha(self.results, self.l, self.l_SM, filepath)
        else:
            writeoutput.results_xml(self.results, self.l, version.__version__, self.dbversion,
                                    filepath)
