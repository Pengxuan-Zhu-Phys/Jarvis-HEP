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

class LilithError(Exception):
    """Base error class for all exceptions raised by Lilith."""

class ExpInputError(LilithError):
    """Exception class in case of error in the reading of the experimental
       input."""

    def __init__(self, filepath, error):
        super(ExpInputError, self).__init__(
            "in file " + filepath + ", " + error)
        self.filepath = filepath

class ExpInputIOError(LilithError):
    """Exception class in case of error when reading the experimental input
       file."""

class UserInputError(LilithError):
    """Exception class in case of error in the reading of the user input."""

class UserInputIOError(UserInputError):
    """Exception class if the user input is read from a file that cannot be
       read."""

class HiggsMassError(UserInputError):
    """Exception class if the Higgs mass is outside the required mass range."""

class LikelihoodComputationError(LilithError):
    """Exception class in case of error when computing the likelihood."""

class ExpNdfComputationError(LilithError):
    """Exception class in case of error when computing the experimental number
       of degrees of freedom."""

class UserMuTotComputationError(LilithError):
    """Exception class in case of error when computing the total signal
       strengths from the individual ones."""

class ReducedCouplingComputationError(LilithError):
    """Exception class in case of error when computing the missing reduced
       couplings."""

class ComputeMuFromReducedCouplingsError(LilithError):
    """Exception class in case of error when computing signal strengths from the
       reduced couplings."""

class OutputError(LilithError):
    """Exception class in case of error when dealing with the output of info."""

class OuputIOError(OutputError):
    """Expection class in case of error when writing any output in a file."""

