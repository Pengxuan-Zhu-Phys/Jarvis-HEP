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

from ..errors import LikelihoodComputationError
import numpy as np

def compute_likelihood(exp_mu, user_mu, user_mode):
    """Computes the likelihood from experimental mu and user mu."""
    likelihood_results = []
    l = 0. # actually -2log(likelihood)
    for mu in exp_mu:
        # compute user mu value scaled to efficiencies
        user_mu_effscaled = {}
        try:
            if mu["dim"] == 1:
                user_mu_effscaled["x"] = 0.
                for (prod,decay),eff_prod in list(mu["eff"]["x"].items()):
                    if mu["sqrts"] not in ["1.96","7","8","7.","8.","7.0","8.0","7+8"]\
                      and (prod == "ggH" or prod == "VBF" or prod == "tHq" or prod == "tHW" or prod == "ggZH"):
                        if user_mode == "reducedcouplings":
                            prod = prod + "13"
                    user_mu_effscaled["x"] += eff_prod*user_mu[prod,decay]
            elif mu["dim"] == 2:
                user_mu_effscaled["x"] = 0.
                for (prod,decay),eff_prod in list(mu["eff"]["x"].items()):
                    if mu["sqrts"] not in ["1.96","7","8","7.","8.","7.0","8.0","7+8"]\
                      and (prod == "ggH" or prod == "VBF" or prod == "tHq" or prod == "tHW" or prod == "ggZH"):
                        if user_mode == "reducedcouplings":
                            prod = prod + "13"
                    user_mu_effscaled["x"] += eff_prod*user_mu[prod,decay]

                user_mu_effscaled["y"] = 0.
                for (prod,decay),eff_prod in list(mu["eff"]["y"].items()):
                    if mu["sqrts"] not in ["1.96","7","8","7.","7.0","8.0","8.","7+8"]\
                      and (prod == "ggH" or prod == "VBF" or prod == "tHq" or prod == "tHW" or prod == "ggZH"):
                        if user_mode == "reducedcouplings":
                            prod = prod + "13"
                    user_mu_effscaled["y"] += eff_prod*user_mu[prod,decay]
            elif mu["dim"] >= 3:
                for i in range(1,mu["dim"]+1):
                    d = "d" + str(i)
                    user_mu_effscaled[d] = 0.
                    for (prod,decay),eff_prod in list(mu["eff"][d].items()):
                        if mu["sqrts"] not in ["1.96","7","8","7.","8.","7.0","8.0","7+8"]\
                          and (prod == "ggH" or prod == "VBF" or prod == "tHq" or prod == "tHW" or prod == "ggZH"):
                            if user_mode == "reducedcouplings":
                                prod = prod + "13"
                        user_mu_effscaled[d] += eff_prod*user_mu[prod,decay]

        except KeyError as s:
            if "s" in ["eff", "x", "y"]:
                # the experimental mu dictionnary is not filled correctly
                raise LikelihoodComputationError(
                    'there are missing elements in exp_mu: key "' + str(s) +
                    '" is not found')
            else:
                # the total user mu dictionnary is not filled correctly
                raise LikelihoodComputationError(
                    'there are missing elements in user_mu_tot: key "' +
                    str(s) + '" is not found')

        try:
            # likelihood computation in case of a type="normal" (odinary Gaussian approximation)
            if mu["type"] == "n":
                if mu["dim"] == 1:
                    if user_mu_effscaled["x"] < mu["bestfit"]["x"]:
                        unc = mu["param"]["uncertainty"]["left"]
                    else:
                        unc = mu["param"]["uncertainty"]["right"]
                    cur_l = ((mu["bestfit"]["x"] - user_mu_effscaled["x"])**2/unc**2)

                elif mu["dim"] == 2:
                    a = mu["param"]["a"]
                    b = mu["param"]["b"]
                    c = mu["param"]["c"]

                    cur_l = a*(mu["bestfit"]["x"] - user_mu_effscaled["x"])**2
                    cur_l += c*(mu["bestfit"]["y"] - user_mu_effscaled["y"])**2
                    cur_l += (2*b*(mu["bestfit"]["x"] - user_mu_effscaled["x"])
                             * (mu["bestfit"]["y"] - user_mu_effscaled["y"]))

                elif mu["dim"] >= 3:
                    mu_vec = np.array([user_mu_effscaled["d1"] - mu["bestfit"]["d1"],
                                       user_mu_effscaled["d2"] - mu["bestfit"]["d2"],
                                       user_mu_effscaled["d3"] - mu["bestfit"]["d3"]])
                    for i in range(4,mu["dim"]+1):
                        d = "d"+str(i)
                        mu_vec = np.append(mu_vec,[user_mu_effscaled[d] - mu["bestfit"][d]])

                    cur_l = mu["param"]["inv_cov_m"].dot(mu_vec).dot(mu_vec.T)

            # likelihood computation in case of a type="variable normal"
            # following "Variable Gaussian 2", Barlow arXiv:physics/0406120v1, Eq. 18
            if mu["type"] == "vn":
                if mu["dim"] == 1:
                    unc_left = abs(mu["param"]["uncertainty"]["left"])
                    unc_right = mu["param"]["uncertainty"]["right"]
                    if unc_left == 0:
                        cur_l = (user_mu_effscaled["x"]-mu["bestfit"]["x"])/unc_right
                    elif unc_right == 0:
                        cur_l = -(user_mu_effscaled["x"]-mu["bestfit"]["x"])/unc_left
                    else:
                        num = user_mu_effscaled["x"] - mu["bestfit"]["x"]
                        den = unc_left*unc_right + (unc_right - unc_left)*num
                        if den == 0:
                            raise LikelihoodComputationError(
                              'divided by zero in 1D Variable Gaussian case')
                        else:
                            cur_l = num**2/den
                if mu["dim"] == 2:
                    p = mu["param"]["correlation"]
                    sig1p = mu["param"]["uncertainty"]["x"]["right"]
                    sig1m = abs(mu["param"]["uncertainty"]["x"]["left"])
                    sig2p = mu["param"]["uncertainty"]["y"]["right"]
                    sig2m = abs(mu["param"]["uncertainty"]["y"]["left"])
                    z10 = mu["bestfit"]["x"]
                    z20 = mu["bestfit"]["y"]
                    z1 = user_mu_effscaled["x"]
                    z2 = user_mu_effscaled["y"]
                    V1 = sig1p*sig1m
                    V1e = sig1p - sig1m
                    V2 = sig2p*sig2m
                    V2e = sig2p - sig2m
                    V1f = V1 + V1e*(z1-z10)
                    V2f = V2 + V2e*(z2-z20)
                    cur_l = 1.0/(1-p**2)*((z1-z10)**2/V1f-2*p*(z1-z10)*(z2-z20)/np.sqrt(V1f*V2f)+(z2-z20)**2/V2f)
                elif mu["dim"] >= 3:
                    mu_vec = np.array([user_mu_effscaled["d1"] - mu["bestfit"]["d1"],
                                       user_mu_effscaled["d2"] - mu["bestfit"]["d2"],
                                       user_mu_effscaled["d3"] - mu["bestfit"]["d3"]])
                    for i in range(4,mu["dim"]+1):
                        d = "d"+str(i)
                        mu_vec = np.append(mu_vec,[user_mu_effscaled[d] - mu["bestfit"][d]])

                    unc_sym = np.sqrt(mu["param"]["VGau"] + mu["param"]["VGau_prime"]*mu_vec)
                    cov_m = unc_sym*mu["param"]["corr_m"]*unc_sym.T
                    inv_cov_m = np.linalg.inv(cov_m)
                    cur_l = inv_cov_m.dot(mu_vec).dot(mu_vec.T)

            # likelihood computation in case of a type="Poisson"
            # following "Generalised Poisson" of Barlow, arXiv:physics/0406120v1, Eq. 10a
            if mu["type"] == "p":
                if mu["dim"] == 1:
                    sigm = abs(mu["param"]["uncertainty"]["left"])
                    sigp = mu["param"]["uncertainty"]["right"]
                    x0 = mu["bestfit"]["x"]
                    x = user_mu_effscaled["x"]

                    gamma = mu["param"]["gamma"]
                    nu = mu["param"]["nu"]
                    alpha = nu*gamma
                    cur_l = -alpha*(x-x0) + nu*np.log(1+alpha*(x-x0)/nu)
                    cur_l = -2.*cur_l

                if mu["dim"] == 2:
                    p = mu["param"]["correlation"]
                    sig1p = mu["param"]["uncertainty"]["x"]["right"]
                    sig1m = abs(mu["param"]["uncertainty"]["x"]["left"])
                    sig2p = mu["param"]["uncertainty"]["y"]["right"]
                    sig2m = abs(mu["param"]["uncertainty"]["y"]["left"])
                    z10 = mu["bestfit"]["x"]
                    z20 = mu["bestfit"]["y"]
                    z1 = user_mu_effscaled["x"]
                    z2 = user_mu_effscaled["y"]

                    gamma1 = mu["param"]["gamma"]["x"]
                    gamma2 = mu["param"]["gamma"]["y"]

                    nu1 = mu["param"]["nu"]["x"]
                    alpha1 = nu1*gamma1
                    nu2 = mu["param"]["nu"]["y"]
                    alpha2 = nu2*gamma2
                    A = mu["param"]["A_corr"]
                    alpha = mu["param"]["alpha_corr"]
                    L2t1 = -alpha1*(z1-z10) + nu1*np.log(1+alpha1*(z1-z10)/nu1)
                    L2t2a = -alpha2*(z2 - z20 + 1/gamma2)*np.exp(alpha*nu1 - A*alpha1*(z1 - z10 + 1/gamma1))
                    L2t2b = -alpha2*(1/gamma2)*np.exp(alpha*nu1 - A*alpha1/gamma1)
                    L2t2c = nu2*np.log(L2t2a/L2t2b)
                    L2t2 = L2t2a - L2t2b + L2t2c
                    cur_l = -2.0*(L2t1 + L2t2)

            # likelihood computation in case of a type="full" (exact likelihood provided in terms of a grid file)
            if mu["type"] == "f":
                if mu["dim"] == 1:
#                    cur_l = mu["Lxy"](user_mu_effscaled["x"]) - mu["LChi2min"]
                    cur_l = mu["Lxy"](user_mu_effscaled["x"])
                elif mu["dim"] == 2:
#                    cur_l = mu["Lxy"](user_mu_effscaled["x"],user_mu_effscaled["y"])[0][0] - mu["LChi2min"]
                    cur_l = mu["Lxy"](user_mu_effscaled["x"],user_mu_effscaled["y"])[0][0]
                if cur_l < 0:
                    cur_l = 0.
        except KeyError as s:
            raise LikelihoodComputationError(
                'there are missing elements in exp_mu: key "' + s +
                '" is not found')

# LDN added a control on the value of loglikelihood
#        if cur_l < 0:
##            if mu["dim"] == 1 and mu["type"] == "f":
##                print("mu =",user_mu_effscaled["x"])
##                print("cur_l min =",mu["LChi2min"])
##                print("cur_l x   =",mu["Lxy"](user_mu_effscaled["x"]))
##                print("cur_l =",cur_l)
#            raise LikelihoodComputationError(
#                'loglikelihood is negative, check the value of mu')
# end


        l += cur_l
        likelihood_results.append(
            {"experiment": mu["experiment"], "source": mu["source"],
            "sqrts": mu["sqrts"], "dim": mu["dim"],
            "type": mu["type"], "eff": mu["eff"], "l": cur_l})

    return likelihood_results, l
