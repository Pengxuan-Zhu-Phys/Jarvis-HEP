// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef GM2_FFUNCTIONS_HPP
#define GM2_FFUNCTIONS_HPP

namespace gm2calc {

/// \f$F_1^C(x)\f$, Eq (54) arXiv:hep-ph/0609168
double F1C(double) noexcept;
/// \f$F_2^C(x)\f$, Eq (55) arXiv:hep-ph/0609168
double F2C(double) noexcept;
/// \f$F_3^C(x)\f$, Eq (37) arXiv:1003.5820
double F3C(double) noexcept;
/// \f$F_4^C(x)\f$, Eq (38) arXiv:1003.5820
double F4C(double) noexcept;
/// \f$F_1^N(x)\f$, Eq (52) arXiv:hep-ph/0609168
double F1N(double) noexcept;
/// \f$F_2^N(x)\f$, Eq (53) arXiv:hep-ph/0609168
double F2N(double) noexcept;
/// \f$F_3^N(x)\f$, Eq (39) arXiv:1003.5820
double F3N(double) noexcept;
/// \f$F_4^N(x)\f$, Eq (40) arXiv:1003.5820
double F4N(double) noexcept;
/// \f$F_a(x)\f$, Eq (6.3a) arXiv:1311.1775
double Fa(double, double) noexcept;
/// \f$F_b(x)\f$, Eq (6.3b) arXiv:1311.1775
double Fb(double, double) noexcept;
/// \f$G_3(x)\f$, Eq (6.4a) arXiv:1311.1775
double G3(double) noexcept;
/// \f$G_4(x)\f$, Eq (6.4b) arXiv:1311.1775
double G4(double) noexcept;
/// \f$I_{abc}(a,b,c)\f$ (arguments are interpreted as unsquared)
double Iabc(double, double, double) noexcept;
/// \f$f_{PS}(z)\f$, Eq (70) arXiv:hep-ph/0609168
double f_PS(double) noexcept;
/// \f$f_S(z)\f$, Eq (71) arXiv:hep-ph/0609168
double f_S(double) noexcept;
/// \f$f_{\tilde{f}}(z)\f$, Eq (72) arXiv:hep-ph/0609168
double f_sferm(double) noexcept;
/// \f$f_l^{H^\pm}(z)\f$, Eq (60) arxiv:1607.06292
double f_CSl(double) noexcept;
/// \f$\mathcal{F}_d^{H^\pm}(x,y,q_u,q_d)\f$, Eq (61) arxiv:1607.06292
double f_CSd(double, double, double, double) noexcept;
/// \f$\mathcal{F}_u^{H^\pm}(x,y,q_u,q_d)\f$, Eq (62) arxiv:1607.06292
double f_CSu(double, double, double, double) noexcept;
/// \f$\mathcal{F}_1(\omega)\f$, Eq (25) arxiv:1502.04199
double F1(double) noexcept;
/// \f$\tilde{\mathcal{F}}_1(\omega)\f$, Eq (26) arxiv:1502.04199
double F1t(double) noexcept;
/// \f$\mathcal{F}_2(\omega)\f$, Eq (27) arxiv:1502.04199
double F2(double) noexcept;
/// \f$\mathcal{F}_3(\omega)\f$, Eq (28) arxiv:1502.04199
double F3(double) noexcept;
/// \f$\tilde{F}_{FZ}(x,y)\f$
double FPZ(double, double) noexcept;
/// \f$F_{FZ}(x,y)\f$
double FSZ(double, double) noexcept;
/// \f$F_{CW}^l(x,y)\f$
double FCWl(double, double) noexcept;
/// \f$F_{CW}^u(x_u,x_d,y_u,y_d,q_u,q_d)\f$
double FCWu(double, double, double, double, double, double) noexcept;
/// \f$F_{CW}^d(x_u,x_d,y_u,y_d,q_u,q_d)\f$
double FCWd(double, double, double, double, double, double) noexcept;
/// \f$\Phi(x,y,z)\f$ with squared masses, Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
double Phi(double x, double y, double z) noexcept;
/// Källén lambda function \f$\lambda^2(x, y, z)\f$
double lambda_2(double x, double y, double z) noexcept;

} // namespace gm2calc

#endif
