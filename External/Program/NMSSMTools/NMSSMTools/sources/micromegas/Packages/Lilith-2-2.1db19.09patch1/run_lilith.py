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

import sys, os, getopt, pprint
from math import floor
from warnings import warn
import lilith

options_short = "hvtsc:m:r:"
options_long = ["help", "verbose", "timer", "silent", "mu=", "couplings=",
                "results="]

def print_usage():
    print("\nUsage of Lilith:")
    print("----------------")
    print("./run_lilith model_input_xml [experimental_input_list]")
    print("Example: ./run_lilith userinput/example_mu.xml [data/latest.list]\n")

def print_options():
    print("[options]")
    print(" -h or --help : dump this help")
    print(" -v or --verbose : print useful information on the screen")
    print(" -t or --timer : print execution time of various tasks")
    print(" -s or --silent : nothing is printed on the screen")
    print(" -m or --mu= : output signal strengths in XML format")
    print((" -c or --couplings= : return missing couplings instead of the " +
           "likelihood"))
    print(" -r file or --results=file : prints the output in a file\n")

def warning(message):
    warn(message, Warning, stacklevel=3)

# at least one argument, check if asking for help
if len(sys.argv) < 2 or (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
    print_usage()
    print_options()
    sys.exit()

user_input_filepath = sys.argv[1]

# default values
exp_input_filepath = ""
verbose = False
timer = False
silent = False

couplingsmode = False
mumode = False
resultsmode = False
signalstrengthmode = False
writtenmu = False
writtencouplings = False
writtenresults = False
optionwarning = ""
resultsfilepath = ""
mufilepath = ""
couplingsfilepath = ""


# at least 2 args: set exp_input_filepath if not an option
if len(sys.argv) >= 3:
    n = 2
    if sys.argv[2][0] != "-":
        n = 3
        exp_input_filepath = sys.argv[2]

    # check all possible options...
    try:
        opts, args = getopt.getopt(sys.argv[n:], options_short,
                                   options_long)
    except getopt.GetoptError as err:
        print(err)
        print_usage()
        print_options()
        sys.exit()

    # and take into account these options
    for o, a in opts:
        if o in ("-v", "--verbose"):
            verbose = True
            if silent:
                print(o, "is incompatible with silent mode.")
                sys.exit()
        elif o in ("-t", "--timer"):
            timer = True
            if silent:
                print(o, "is incompatible with silent mode.")
                sys.exit()
        elif o in ("-s", "--silent"):
            silent = True
            if verbose or timer:
                print(o, "is incompatible with verbose mode.")
                sys.exit()
        elif o in ("-h", "--help"):
            print_usage()
            print_options()
            sys.exit()
        elif o in ("-m", "--mu"):
            mumode = True
            mufilepath = a
        elif o in ("-c", "--couplings"):
            couplingsmode = True
            couplingsfilepath = a
        elif o in ("-r", "--results"):
            resultsmode = True
            resultsfilepath = a
        else:
            print("unhandled option")
            sys.exit()


# turning to the calculation of the likelihood
Lilithcalc = lilith.Lilith(verbose, timer)


if exp_input_filepath:
    Lilithcalc.readexpinput(exp_input_filepath)
else:
    Lilithcalc.readexpinput()
Lilithcalc.readuserinputfile(user_input_filepath)
Lilithcalc.computelikelihood(userfilepath=user_input_filepath)
Lilithcalc.computecouplings()

if couplingsmode and not Lilithcalc.couplings:
    signalstrengthmode = True
    if not silent:
        optionwarning="Cannot compute couplings in signal strengths mode"\
                    +", -c/--couplings option will be ignored\n"


if not silent:
    lv = len(lilith.__version__)
    n_lspace = 0
    if lv > 0 and lv < 11:
        n_lspace = int(floor(0.5*(11-lv)))
    
    lv2 = len(Lilithcalc.dbversion)
    n_lspace2 = 0
    if lv2 > 0 and lv2 < 11:
        n_lspace2 = int(floor(0.5*(11-lv2)))

    print(("\n" +
           "                           " +
           "<><><><><><><><><><><><><>"))
    print(("                           " +
           " "*n_lspace + "Lilith version " + lilith.__version__+
           "\n                          "+" "*n_lspace2 +
           "database version " + Lilithcalc.dbversion))
    print(("                           " +
           "<><><><><><><><><><><><><>\n"))

    print(". User input: " + user_input_filepath + "\n")

    if optionwarning:
        warning(optionwarning)

    if not couplingsmode:
        if not exp_input_filepath:
            print((". Experimental input: latest LHC results " +
                   "[data/latest.list]"))
        else:
            print(". Experimental input: " + exp_input_filepath)
        print(("                      (" + str(len(Lilithcalc.exp_mu)) +
               " files, Ndof = " + str(Lilithcalc.exp_ndf) + ")"))
        print("-2log(likelihood) = " + str(round(Lilithcalc.l,6)))
    else:
        if optionwarning=="":
            pp = pprint.PrettyPrinter()
            print("Reduced couplings are:")
            if len(Lilithcalc.couplings) == 1:
                pp.pprint(Lilithcalc.couplings[0])
            else:
                pp.pprint(Lilithcalc.couplings)


if couplingsmode and couplingsfilepath != "" and not signalstrengthmode:
    Lilithcalc.writecouplings(couplingsfilepath)
    writtencouplings = True

if mumode and mufilepath !="":
    Lilithcalc.writesignalstrengths(mufilepath)
    writtenmu = True

if resultsmode and resultsfilepath !="":
    if resultsfilepath.split(".")[-1].lower() == "slha":
        Lilithcalc.computeSMlikelihood(userfilepath=user_input_filepath)
        Lilithcalc.writeresults(resultsfilepath, slha=True)
        writtenresults = True
    else:
        Lilithcalc.writeresults(resultsfilepath, slha=False)
        writtenresults = True


if writtenmu and not silent:
        print(". Signal strengths have been written in "+mufilepath)
if writtencouplings and not silent:
        print(". Couplings have been written in "+couplingsfilepath)
if writtenresults and not silent:
        print(". Results have been written in "+resultsfilepath)

