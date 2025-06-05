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

#ifndef GM2_CONSTANTS_HPP
#define GM2_CONSTANTS_HPP

namespace gm2calc {

// fine-structure constant in the Thompson limit (Q = 0) from PDG (2019)
constexpr double ALPHA_EM_THOMPSON = 1.0/137.035999084;

// quark and lepton contributions to the on-shell renormalized photon
// vacuum polarization
constexpr double DELTA_ALPHA_EM_MZ =
      + 0.0314979  // leptonic
      - 0.00007180 // top
      + 0.027611;  // hadronic, arXiv:1802.02995

// fine-structure constant at Q = MZ
constexpr double ALPHA_EM_MZ = ALPHA_EM_THOMPSON / (1 - DELTA_ALPHA_EM_MZ);

// strong coupling alpha_s at Q = MZ
constexpr double ALPHA_S_MZ = 0.1184;

// Higgs boson pole mass
constexpr double MH = 125.09;

// W boson pole mass
constexpr double MW = 80.385;

// Z boson pole mass
constexpr double MZ = 91.1876;

// Up quark mass
constexpr double MU = 0.0022;

// Charm quark mass
constexpr double MC = 1.28;

// Top quark pole mass
constexpr double MT = 173.34;

// Down quark mass
constexpr double MD = 0.0047;

// Strange quark mass
constexpr double MS = 0.096;

// SM MS-bar bottom quark mass mb at the scale mb
constexpr double MBMB = 4.18;

// Electron pole mass
constexpr double ME = 0.000510998928;

// Muon pole mass
constexpr double MM = 0.1056583715;

// Tau lepton pole mass
constexpr double ML = 1.777;

// From Vus/Vud in global CKM fit, PDG
constexpr double CKM_THETA12 = 0.229206;

// From Vub in global CKM fit, PDG
constexpr double CKM_THETA13 = 0.003960;

// From Vcb/Vtb in global CKM fit, PDG
constexpr double CKM_THETA23 = 0.042223;

constexpr double CKM_DELTA = 0;

} // namespace gm2calc

#endif
