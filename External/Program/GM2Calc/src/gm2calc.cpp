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

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_version.h"
#include "gm2calc/MSSMNoFV_onshell.hpp"
#include "gm2calc/THDM.hpp"

#include "MSSMNoFV/gm2_1loop_helpers.hpp"
#include "MSSMNoFV/gm2_2loop_helpers.hpp"
#include "THDM/gm2_2loop_helpers.hpp"
#include "gm2_config_options.hpp"
#include "gm2_log.hpp"
#include "gm2_slha_io.hpp"

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include "json.hpp"
#include <fstream>

using json = nlohmann::json;

#define FORMAT_AMU(amu) std::scientific << std::setprecision(8) << std::setw(15) << (amu)
#define FORMAT_DEL(amu) std::scientific << std::setprecision(8) << std::setw(14) << (amu)
#define FORMAT_PCT(pct) std::fixed << std::setprecision(1) << std::setw(2) << (pct)

namespace {

/// MSSMNoFV reader function
using MSSMNoFV_reader = std::function<void(
   gm2calc::MSSMNoFV_onshell&, const gm2calc::GM2_slha_io& slha_io)>;

/// MSSMNoFV writer function
using MSSMNoFV_writer = std::function<void(
   const gm2calc::MSSMNoFV_onshell&, const gm2calc::Config_options& options,
   gm2calc::GM2_slha_io& slha_io)>;

/// THDM writer function
using THDM_writer = std::function<void(const gm2calc::THDM&,
                                       const gm2calc::Config_options& options,
                                       gm2calc::GM2_slha_io& slha_io)>;

/**
 * @class Gm2_cmd_line_options
 * @brief command line options for GM2Calc
 */
struct Gm2_cmd_line_options {
   enum E_input_type { SLHA, GM2Calc, THDM };

   std::string input_source; ///< input source (file name or `-' for stdin)
   E_input_type input_type{SLHA}; ///< input format (SLHA, GM2Calc or THDM)

   static bool starts_with(const std::string& str, const std::string& prefix) {
      return str.compare(0, prefix.size(), prefix) == 0;
   }
};

void print_usage(const char* program_name)
{
   std::cout <<
      "Usage: " << program_name << " [options]\n"
      "Options:\n"
      "  --slha-input-file=<source>      SLHA input source (file name or - for stdin)\n"
      "  --gm2calc-input-file=<source>   GM2Calc input source (file name or - for stdin)\n"
      "  --thdm-input-file=<source>      THDM input source (file name or - for stdin)\n"
      "  --help,-h                       print this help message\n"
      "  --version,-v                    print version number"
      "\n";
}

/**
 * Checks whether the given string is a specification of an input file
 * and processes it, if so.
 *
 * @param str option string
 * @param options option struct to modify
 *
 * @return true if option string has been processed, false otherwise
 */
bool process_input_type(const std::string& str, Gm2_cmd_line_options& options)
{
   static const struct Input_file_options {
      std::string input_source;
      Gm2_cmd_line_options::E_input_type input_type;
   } input_file_options[] = {
      { "--slha-input-file="   , Gm2_cmd_line_options::SLHA    },
      { "--gm2calc-input-file=", Gm2_cmd_line_options::GM2Calc },
      { "--thdm-input-file="   , Gm2_cmd_line_options::THDM    }
   };

   for (const auto& ifo: input_file_options) {
      if (Gm2_cmd_line_options::starts_with(str, ifo.input_source)) {
         options.input_source = str.substr(ifo.input_source.length());
         options.input_type = ifo.input_type;
         return true; // option has been processed
      }
   }

   return false; // option has not been processed
}

/**
 * Parses command line options
 *
 * @param argc number of command line arguments
 * @param argv array of command line arguments
 *
 * @return object with extracted information
 */
Gm2_cmd_line_options get_cmd_line_options(int argc, const char* argv[])
{
   Gm2_cmd_line_options options;

   for (int i = 1; i < argc; ++i) {
      const std::string option_string(argv[i]);

      if (process_input_type(option_string, options)) {
         continue;
      }

      if (option_string == "--help" || option_string == "-h") {
         print_usage(argv[0]);
         exit(EXIT_SUCCESS);
      }

      if (option_string == "--version" || option_string == "-v") {
         std::cout << GM2CALC_VERSION << '\n';
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option_string);
      exit(EXIT_FAILURE);
   }

   return options;
}

/**
 * Set the config options to default values, depending on the input
 * parameter set (chosen by the user).
 *
 * If SLHA input format has been selected, the default output format
 * will be SLHA format.  By default \f$a_\mu\f$ will be written to the
 * block GM2Calc[0].
 *
 * If GM2Calc input format has been chosen, the default values set in
 * \a Config_options are used.
 *
 * @param config_options configuration options
 * @param options command line options
 */
void set_to_default(gm2calc::Config_options& config_options,
                    const Gm2_cmd_line_options& options)
{
   switch (options.input_type) {
   case Gm2_cmd_line_options::SLHA:
      config_options.output_format = gm2calc::Config_options::GM2Calc;
      break;
   case Gm2_cmd_line_options::GM2Calc:
      config_options.output_format = gm2calc::Config_options::Detailed;
      break;
   case Gm2_cmd_line_options::THDM:
      config_options.output_format = gm2calc::Config_options::GM2Calc;
      break;
   default:
      throw gm2calc::ESetupError("Unknown input option");
      break;
   }
}

/**
 * Prints output if an error has occured.
 *
 * @param error error object
 * @param slha_io SLHA object
 * @param config_options configuration options
 */
void print_error(const gm2calc::Error& error,
                 gm2calc::GM2_slha_io& slha_io,
                 const gm2calc::Config_options& config_options)
{
   switch (config_options.output_format) {
   case gm2calc::Config_options::NMSSMTools:
   case gm2calc::Config_options::SPheno:
   case gm2calc::Config_options::GM2Calc:
      // print SPINFO block with error description
      slha_io.fill_block_entry("SPINFO", 1, "GM2Calc");
      slha_io.fill_block_entry("SPINFO", 2, GM2CALC_VERSION);
      slha_io.fill_block_entry("SPINFO", 4, error.what());
      slha_io.write_to_stream(std::cout);
      break;
   default:
      ERROR(error.what());
      break;
   }
}

/**
 * Calculates a_mu for a given set of configuration options (loop
 * order, tan(beta) resummation).
 *
 * @param model model (must be initialized)
 * @param options configuration options
 *
 * @return a_mu
 */
double calculate_amu(const gm2calc::MSSMNoFV_onshell& model,
                     const gm2calc::Config_options& options)
{
   double result = 0.0;

   if (options.tanb_resummation) {
      if (options.loop_order > 0) {
         result += gm2calc::calculate_amu_1loop(model);
      }
      if (options.loop_order > 1) {
         result += gm2calc::calculate_amu_2loop(model);
      }
   } else {
      // no tan(beta) resummation
      if (options.loop_order > 0) {
         result += gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model);
      }
      if (options.loop_order > 1) {
         result += gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model);
      }
   }

   return result;
}

/**
 * Calculates a_mu for a given set of configuration options (loop
 * order).
 *
 * @param model model (must be initialized)
 * @param options configuration options
 *
 * @return a_mu
 */
double calculate_amu(const gm2calc::THDM& model,
                     const gm2calc::Config_options& options)
{
   double result = 0.0;

   if (options.loop_order > 0) {
      result += gm2calc::calculate_amu_1loop(model);
   }
   if (options.loop_order > 1) {
      result += gm2calc::calculate_amu_2loop(model);
   }

   return result;
}

/**
 * Calculates uncertainty of a_mu for a given set of configuration
 * options (loop order, tan(beta) resummation).
 *
 * @param model model (must be initialized)
 * @param options configuration options
 *
 * @return a_mu
 */
template<class Model>
double calculate_uncertainty(const Model& model,
                             const gm2calc::Config_options& options)
{
   double result = std::numeric_limits<double>::signaling_NaN();

   switch (options.loop_order) {
   case 0:
      result = gm2calc::calculate_uncertainty_amu_0loop(model);
      break;
   case 1:
      result = gm2calc::calculate_uncertainty_amu_1loop(model);
      break;
   case 2:
      result = gm2calc::calculate_uncertainty_amu_2loop(model);
      break;
   default:
      ERROR("loop order > 2 not supported!");
      break;
   }

   return result;
}

/**
 * Reads parameters from SLHA i/o object (SLHA scheme) and initializes
 * model accordingly.
 *
 * @param model model to initialize
 * @param slha_io SLHA i/o object to read parameters from
 */
struct SLHA_reader {
   void operator()(gm2calc::MSSMNoFV_onshell& model,
                   const gm2calc::GM2_slha_io& slha_io)
   {
      slha_io.fill_slha(model);
      model.convert_to_onshell();
   }
};

/**
 * Reads parameters from SLHA i/o object in GM2Calc-specific input
 * scheme and initializes model accordingly.
 *
 * @param model model to initialize
 * @param slha_io SLHA i/o object to read parameters from
 */

struct GM2Calc_reader {
   void operator()(gm2calc::MSSMNoFV_onshell& model,
                   const gm2calc::GM2_slha_io& slha_io)
   {
      slha_io.fill_gm2calc(model);
      model.calculate_masses();
   }
};

/**
 * Reads THDM parameters from SLHA i/o object and initializes model
 * accordingly.
 *
 * @param model model to initialize
 * @param slha_io SLHA i/o object to read parameters from
 */

struct THDM_reader {
   gm2calc::THDM operator()(const gm2calc::GM2_slha_io& slha_io,
                            const gm2calc::Config_options& options)
   {
      gm2calc::SM sm;
      gm2calc::thdm::Mass_basis mass_basis;
      gm2calc::thdm::Gauge_basis gauge_basis;
      slha_io.fill(sm);
      slha_io.fill(mass_basis);
      slha_io.fill(gauge_basis);

      gm2calc::thdm::Config thdm_config;
      thdm_config.force_output = options.force_output;
      thdm_config.running_couplings = options.running_couplings;

      // test for unset parameters to decide which basis to use
      if ((mass_basis.mh != 0 || mass_basis.mH != 0 ||
           mass_basis.mA != 0 || mass_basis.mHp != 0 ||
           mass_basis.sin_beta_minus_alpha != 0) &&
          gauge_basis.lambda.head<5>().cwiseAbs().maxCoeff() == 0) {
         return gm2calc::THDM(mass_basis, sm, thdm_config);
      } else if (mass_basis.mh == 0 && mass_basis.mH == 0 &&
                 mass_basis.mA == 0 && mass_basis.mHp == 0 &&
                 mass_basis.sin_beta_minus_alpha == 0 &&
                 gauge_basis.lambda.head<5>().cwiseAbs().maxCoeff() != 0) {
         return gm2calc::THDM(gauge_basis, sm, thdm_config);
      } else {
         throw gm2calc::EInvalidInput("Cannot distinguish between mass and gauge basis.");
      }
   }
};

/**
 * Prints a_mu (or the uncertainty) to stdout.
 *
 * @param model the model (must be initialized)
 * @param options calculation options
 * @param slha_io SLHA i/o object where results are stored
 */
template<class Model>
struct Minimal_writer {
   void operator()(const Model& model,
                   const gm2calc::Config_options& options,
                   gm2calc::GM2_slha_io& /* unused */)
   {
      const auto value = options.calculate_uncertainty
                            ? calculate_uncertainty(model, options)
                            : calculate_amu(model, options);

      std::cout << std::scientific << std::setprecision(8) << value << '\n';

      // Write JSON output using nlohmann::json
      json root;
      root["a_mu"] = value;
      std::ofstream ofs("gm2calc_amu.json");
      if (ofs.is_open()) {
         ofs << std::setw(2) << root << std::endl;
         std::cout << "gm2calc_amu.json generated\n";
      } else {
         std::cerr << "Error: cannot open gm2calc_amu.json for writing\n";
      }
   }
};

template<class Model>
struct Detailed_writer;

/**
 * Prints detailed a_mu calculation (1-loop w/ and w/o tan(beta)
 * resummation, 2-loop, and different contributions).
 *
 * @param model the model (must be initialized)
 * @param options calculation options
 * @param slha_io SLHA i/o object where results are stored
 */
template<>
struct Detailed_writer<gm2calc::MSSMNoFV_onshell> {
   void operator()(const gm2calc::MSSMNoFV_onshell& model,
                   const gm2calc::Config_options& /* unused */,
                   gm2calc::GM2_slha_io& /* unused */)
   {
      const std::string error_str = model.get_problems().have_problem()
                                       ? model.get_problems().get_problems() +
                                            " (with tan(beta) resummation)\n\n"
                                       : "";

      const double amu_1l = gm2calc::calculate_amu_1loop(model);
      const double amu_2l_photonic_chipm = gm2calc::amu2LChipmPhotonic(model);
      const double amu_2l_photonic_chi0 = gm2calc::amu2LChi0Photonic(model);
      const double amu_2l_a_sfermion = gm2calc::amu2LaSferm(model);
      const double amu_2l_a_cha = gm2calc::amu2LaCha(model);
      const double amu_2l_ferm_sferm_approx = gm2calc::amu2LFSfapprox(model);
      const double amu_2l = gm2calc::calculate_amu_2loop(model);
      const double amu_2l_uncertainty =
         gm2calc::calculate_uncertainty_amu_2loop(model);
      const double tan_beta_cor = gm2calc::tan_beta_cor(model);

      // no tan(beta) resummation
      double amu_1l_non_tan_beta_resummed = 0.;
      double amu_2l_non_tan_beta_resummed = 0.;
      std::string error_str_non_tan_beta_resummation;

      try {
         // w/o tan(beta) resummation, allow throwing exceptions
         gm2calc::MSSMNoFV_onshell model_except(model);
         model_except.do_force_output(false);
         amu_1l_non_tan_beta_resummed =
            gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model_except);
         amu_2l_non_tan_beta_resummed =
            gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model_except);
      } catch (const gm2calc::Error& error) {
         error_str_non_tan_beta_resummation =
            " (" + std::string(error.what()) + ")";
         // try to redo calculation w/o throwing an exception
         gm2calc::MSSMNoFV_onshell model_no_except(model);
         model_no_except.do_force_output(true);
         amu_1l_non_tan_beta_resummed =
            gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model_no_except);
         amu_2l_non_tan_beta_resummed =
            gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model_no_except);
      }

      const double amu_2l_tanb_approx =
         (tan_beta_cor - 1.) * amu_1l_non_tan_beta_resummed;

      const double amu_best = amu_1l + amu_2l;

      std::cout
         << "====================================================================\n"
            "   amu (1-loop + 2-loop best) = "
         << FORMAT_AMU(amu_best) << " +- "
         << FORMAT_DEL(amu_2l_uncertainty) << '\n'
         << "====================================================================\n"
            "\n"
         << error_str
         << "==============================\n"
            "   amu (1-loop) corrections\n"
            "==============================\n"
            "\n"
            "full 1L with tan(beta) resummation:\n"
            "   chi^0     " << FORMAT_AMU(gm2calc::amu1LChi0(model)) << '\n'
         << "   chi^+-    " << FORMAT_AMU(gm2calc::amu1LChipm(model)) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(amu_1l)
                            << " (" << FORMAT_PCT(100. * amu_1l / amu_best)
                            << "% of full 1L + 2L result)\n"
            "\n"
            "full 1L without tan(beta) resummation:\n"
            "             " << FORMAT_AMU(amu_1l_non_tan_beta_resummed)
         << error_str_non_tan_beta_resummation << '\n'
         << "\n"
            "1L approximation with tan(beta) resummation:\n"
            "   W-H-nu    " << FORMAT_AMU(gm2calc::amu1LWHnu(model) * tan_beta_cor) << '\n'
         << "   W-H-muL   " << FORMAT_AMU(gm2calc::amu1LWHmuL(model) * tan_beta_cor) << '\n'
         << "   B-H-muL   " << FORMAT_AMU(gm2calc::amu1LBHmuL(model) * tan_beta_cor) << '\n'
         << "   B-H-muR   " << FORMAT_AMU(gm2calc::amu1LBHmuR(model) * tan_beta_cor) << '\n'
         << "   B-muL-muR " << FORMAT_AMU(gm2calc::amu1LBmuLmuR(model) * tan_beta_cor) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(gm2calc::amu1Lapprox(model)) << '\n'
         << "\n"
            "==============================\n"
            "   amu (2-loop) corrections\n"
            "==============================\n"
            "\n"
            "2L best with tan(beta) resummation:\n"
            "             " << FORMAT_AMU(amu_2l) << " (" << FORMAT_PCT(100. * amu_2l / amu_best)
         << "% of full 1L + 2L result)\n"
            "\n"
            "2L best without tan(beta) resummation:\n"
            "             " << FORMAT_AMU(amu_2l_non_tan_beta_resummed)
         << error_str_non_tan_beta_resummation << '\n'
         << "\n"
            "photonic with tan(beta) resummation:\n"
            "   chi^0     " << FORMAT_AMU(amu_2l_photonic_chi0) << '\n'
         << "   chi^+-    " << FORMAT_AMU(amu_2l_photonic_chipm) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(amu_2l_photonic_chipm + amu_2l_photonic_chi0)
                            << " (" << FORMAT_PCT(100. * (amu_2l_photonic_chipm + amu_2l_photonic_chi0) / amu_best)
                            << "% of full 1L + 2L result)\n"
            "\n"
            "fermion/sfermion approximation with tan(beta) resummation:\n"
            "   W-H-nu    " << FORMAT_AMU(gm2calc::amu2LWHnu(model) * tan_beta_cor) << '\n'
         << "   W-H-muL   " << FORMAT_AMU(gm2calc::amu2LWHmuL(model) * tan_beta_cor) << '\n'
         << "   B-H-muL   " << FORMAT_AMU(gm2calc::amu2LBHmuL(model) * tan_beta_cor) << '\n'
         << "   B-H-muR   " << FORMAT_AMU(gm2calc::amu2LBHmuR(model) * tan_beta_cor) << '\n'
         << "   B-muL-muR " << FORMAT_AMU(gm2calc::amu2LBmuLmuR(model) * tan_beta_cor) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(amu_2l_ferm_sferm_approx)
                            << " (" << FORMAT_PCT(100. * amu_2l_ferm_sferm_approx / amu_best)
                            << "% of full 1L + 2L result)\n"
            "\n"
            "2L(a) (1L insertions into 1L SM diagram) with tan(beta) "
            "resummation:\n"
            "   sfermion  " << FORMAT_AMU(amu_2l_a_sfermion) << '\n'
         << "   cha^+-    " << FORMAT_AMU(amu_2l_a_cha) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(amu_2l_a_sfermion + amu_2l_a_cha)
                            << " (" << FORMAT_PCT(100. * (amu_2l_a_sfermion + amu_2l_a_cha) / amu_best)
                            << "% of full 1L + 2L result)\n"
            "\n"
            "tan(beta) correction:\n"
            "   amu(1L) * (1 / (1 + Delta_mu) - 1) = " << FORMAT_AMU(amu_2l_tanb_approx)
                            << " (" << FORMAT_PCT(100. * amu_2l_tanb_approx / amu_1l_non_tan_beta_resummed)
         << "%)\n";

      // Write detailed JSON output using nlohmann::json
      {
          json root;
          root["amu_1l"] = amu_1l;
          root["amu_2l_photonic_chipm"] = amu_2l_photonic_chipm;
          root["amu_2l_photonic_chi0"] = amu_2l_photonic_chi0;
          root["amu_2l_a_sfermion"] = amu_2l_a_sfermion;
          root["amu_2l_a_cha"] = amu_2l_a_cha;
          root["amu_2l_ferm_sferm_approx"] = amu_2l_ferm_sferm_approx;
          root["amu_2l"] = amu_2l;
          root["amu_2l_uncertainty"] = amu_2l_uncertainty;
          root["amu_2l_tanb_approx"] = amu_2l_tanb_approx;
          root["amu_1l_non_tan_beta_resummed"] = amu_1l_non_tan_beta_resummed;
          root["amu_2l_non_tan_beta_resummed"] = amu_2l_non_tan_beta_resummed;
          root["amu_best"] = amu_best;
          root["tan_beta_correction"] = tan_beta_cor;
          if (!error_str.empty()) {
              root["error_str"] = error_str;
          }
          if (!error_str_non_tan_beta_resummation.empty()) {
              root["error_str_non_tan_beta_resummation"] = error_str_non_tan_beta_resummation;
          }
          std::ofstream ofs_det("gm2calc_detailed.json");
          if (ofs_det.is_open()) {
              ofs_det << std::setw(2) << root << std::endl;
              std::cout << "gm2calc_detailed.json generated\n";
          } else {
              std::cerr << "Error: cannot open gm2calc_detailed.json for writing\n";
          }
      }
   }
};

/**
 * Prints detailed a_mu calculation (1-loop, 2-loop, and different
 * contributions).
 *
 * @param model the model (must be initialized)
 * @param options calculation options
 * @param slha_io SLHA i/o object where results are stored
 */
template<>
struct Detailed_writer<gm2calc::THDM> {
   void operator()(const gm2calc::THDM& model,
                   const gm2calc::Config_options& /* unused */,
                   gm2calc::GM2_slha_io& /* unused */)
   {
      const double amu_1l = gm2calc::calculate_amu_1loop(model);
      const double amu_2l = gm2calc::calculate_amu_2loop(model);
      const double amu_2l_B = gm2calc::calculate_amu_2loop_bosonic(model);
      const double amu_2l_F = gm2calc::calculate_amu_2loop_fermionic(model);
      const double amu_2l_uncertainty = gm2calc::calculate_uncertainty_amu_2loop(model);
      const double amu_best = amu_1l + amu_2l;

      std::cout
         << "====================================================================\n"
            "   amu (1-loop + 2-loop) = "
         << FORMAT_AMU(amu_best) << " +- "
         << FORMAT_DEL(amu_2l_uncertainty) << '\n'
         << "====================================================================\n"
            "\n"
         << "==============================\n"
            "   amu (1-loop) corrections\n"
            "==============================\n"
            "\n"
            "full 1L: " << FORMAT_AMU(amu_1l)
                        << " (" << FORMAT_PCT(100. * amu_1l / amu_best)
                        << "% of full 1L + 2L result)\n"
            "\n"
            "==============================\n"
            "   amu (2-loop) corrections\n"
            "==============================\n"
            "\n"
            "bosonic   2L: " << FORMAT_AMU(amu_2l_B) << " (" << FORMAT_PCT(100. * amu_2l_B / amu_2l) << "% of 2L result)\n"
         << "fermionic 2L: " << FORMAT_AMU(amu_2l_F) << " (" << FORMAT_PCT(100. * amu_2l_F / amu_2l) << "% of 2L result)\n"
         << "sum         : " << FORMAT_AMU(amu_2l) << " (" << FORMAT_PCT(100. * amu_2l / amu_best)
         << "% of full 1L + 2L result)\n";

      // Write detailed JSON output for THDM using nlohmann::json
      {
          json root;
          root["amu_1l"] = amu_1l;
          root["amu_2l"] = amu_2l;
          root["amu_2l_bosonic"] = amu_2l_B;
          root["amu_2l_fermionic"] = amu_2l_F;
          root["amu_2l_uncertainty"] = amu_2l_uncertainty;
          root["amu_best"] = amu_best;
          std::ofstream ofs_det("gm2calc_detailed_thdm.json");
          if (ofs_det.is_open()) {
              ofs_det << std::setw(2) << root << std::endl;
              std::cout << "gm2calc_detailed_thdm.json generated\n";
          } else {
              std::cerr << "Error: cannot open gm2calc_detailed_thdm.json for writing\n";
          }
      }
   }
};

/**
 * Calculates a_mu (and potentially also the uncertainty) and writes
 * it to the SLHA i/o object.
 *
 * @param model the model (must be initialized)
 * @param options calculation options
 * @param slha_io SLHA i/o object where results are stored
 */
template<class Model>
struct SLHA_writer {
   /// SLHA entry (block name, key, value, comment)
   using SLHA_entry = std::tuple<std::string, int, double, std::string>;

   void operator()(const Model& model,
                   const gm2calc::Config_options& options,
                   gm2calc::GM2_slha_io& slha_io)
   {
      const SLHA_entry amu_entry = [&] {
         const auto amu = calculate_amu(model, options);
         const auto amu_comment = "Delta(g-2)_muon/2";

         switch (options.output_format) {
         case gm2calc::Config_options::NMSSMTools:
            return SLHA_entry{"LOWEN", 6, amu, amu_comment};
         case gm2calc::Config_options::SPheno:
            return SLHA_entry{"SPhenoLowEnergy", 21, amu, amu_comment};
         default:
            break;
         }

         return SLHA_entry{"GM2CalcOutput", 0, amu, amu_comment};
      }();

      set_SLHA_value(slha_io, amu_entry);

      if (options.calculate_uncertainty) {
         const auto damu = calculate_uncertainty(model, options);
         const SLHA_entry damu_entry{"GM2CalcOutput", 1, damu,
                                     "uncertainty of Delta(g-2)_muon/2"};
         set_SLHA_value(slha_io, damu_entry);
      }

      if (model.get_problems().have_warning()) {
         slha_io.fill_block_entry("SPINFO", 1, "GM2Calc");
         slha_io.fill_block_entry("SPINFO", 2, GM2CALC_VERSION);
         slha_io.fill_block_entry("SPINFO", 3, model.get_problems().get_warnings());
      }

      slha_io.write_to_stream(std::cout);
   }

private:
   /**
    * Sets a entry in a given SLHA block/key.
    *
    * @param slha_io SLHA input/output object
    * @param entry tuple defining the SLHA (block name, key, value, comment)
    */
   void set_SLHA_value(gm2calc::GM2_slha_io& slha_io, const SLHA_entry& entry)
   {
      slha_io.fill_block_entry(std::get<0>(entry), std::get<1>(entry),
                               std::get<2>(entry), std::get<3>(entry));
   }
};

/**
 * Class which handles input/output for the MSSM.
 */
class MSSMNoFV_setup
{
public:
   MSSMNoFV_setup(const gm2calc::Config_options& options_,
                  const MSSMNoFV_reader& reader_,
                  const MSSMNoFV_writer& writer_)
      : options(options_), reader(reader_), writer(writer_)
   {
   }

   /// read from SLHA and write to output
   int run(gm2calc::GM2_slha_io& slha_io)
   {
      if (!reader) {
         throw gm2calc::ESetupError("No reader set");
      }
      if (!writer) {
         throw gm2calc::ESetupError("No writer set");
      }

      gm2calc::MSSMNoFV_onshell model;
      model.do_force_output(options.force_output);
      model.set_verbose_output(options.verbose_output);

      reader(model, slha_io);

      if (options.verbose_output) {
         VERBOSE(model);
      }

      if (model.get_problems().have_problem() ||
          model.get_problems().have_warning()) {
         std::cerr << model.get_problems() << '\n';
      }

      writer(model, options, slha_io);

      return model.get_problems().have_problem() ? EXIT_FAILURE : EXIT_SUCCESS;
   }

private:
   gm2calc::Config_options options;
   MSSMNoFV_reader reader{nullptr};
   MSSMNoFV_writer writer{nullptr};
};

/**
 * Class which handles input/output for the MSSM.
 */
class THDM_setup
{
public:
   THDM_setup(const gm2calc::Config_options& options_,
              const THDM_reader& reader_,
              const THDM_writer& writer_)
      : options(options_), reader(reader_), writer(writer_)
   {
   }

   /// read from SLHA and write to output
   int run(gm2calc::GM2_slha_io& slha_io)
   {
      if (!writer) {
         throw gm2calc::ESetupError("No writer set");
      }

      gm2calc::THDM model = reader(slha_io, options);

      if (options.verbose_output) {
         VERBOSE(model);
      }

      writer(model, options, slha_io);

      return EXIT_SUCCESS;
   }

private:
   gm2calc::Config_options options;
   THDM_reader reader;
   THDM_writer writer;
};

/**
 * Returns properly configured (but not initialized) MSSMNoFV_setup object.
 *
 * @param input_type type of input (SLHA/GM2Calc)
 * @param options configuration options
 *
 * @return MSSMNoFV_setup object
 */
MSSMNoFV_setup make_mssmnofv_setup(
   Gm2_cmd_line_options::E_input_type input_type,
   const gm2calc::Config_options& options)
{
   const auto reader = [&] () -> MSSMNoFV_reader {
      switch (input_type) {
      case Gm2_cmd_line_options::SLHA:
         return SLHA_reader();
      case Gm2_cmd_line_options::GM2Calc:
         return GM2Calc_reader();
      case Gm2_cmd_line_options::THDM:
         break;
      }
      throw gm2calc::ESetupError("Unknown input type");
   }();

   const auto writer = [&] () -> MSSMNoFV_writer {
      switch (options.output_format) {
      case gm2calc::Config_options::Minimal:
         return Minimal_writer<gm2calc::MSSMNoFV_onshell>();
      case gm2calc::Config_options::Detailed:
         return Detailed_writer<gm2calc::MSSMNoFV_onshell>();
      default:
         break;
      }
      return SLHA_writer<gm2calc::MSSMNoFV_onshell>();
   }();

   return MSSMNoFV_setup(options, reader, writer);
}

/**
 * Returns properly configured (but not initialized) MSSMNoFV_setup object.
 *
 * @param input_type type of input (SLHA/GM2Calc)
 * @param options configuration options
 *
 * @return MSSMNoFV_setup object
 */
THDM_setup make_thdm_setup(const gm2calc::Config_options& options)
{
   const auto writer = [&] () -> THDM_writer {
      switch (options.output_format) {
      case gm2calc::Config_options::Minimal:
         return Minimal_writer<gm2calc::THDM>();
      case gm2calc::Config_options::Detailed:
         return Detailed_writer<gm2calc::THDM>();
      default:
         break;
      }
      return SLHA_writer<gm2calc::THDM>();
   }();

   return THDM_setup(options, THDM_reader(), writer);
}

} // anonymous namespace

int main(int argc, const char* argv[])
{
   Gm2_cmd_line_options options(get_cmd_line_options(argc, argv));

   if (options.input_source.empty()) {
      ERROR(std::string("No input source given!\n") +
            "Examples: \n" +
            "   " + argv[0] + " --slha-input-file=<file>     # MSSM SLHA input\n"
            "   " + argv[0] + " --gm2calc-input-file=<file>  # MSSM GM2Calc input\n"
            "   " + argv[0] + " --thdm-input-file=<file>     # THDM input");
      return EXIT_FAILURE;
   }

   gm2calc::GM2_slha_io slha_io;
   gm2calc::Config_options config_options;
   int exit_code = EXIT_SUCCESS;

   try {
      set_to_default(config_options, options);
      slha_io.read_from_source(options.input_source);
      slha_io.fill(config_options);

      switch (options.input_type) {
      case Gm2_cmd_line_options::SLHA:
      case Gm2_cmd_line_options::GM2Calc: {
         auto setup = make_mssmnofv_setup(options.input_type, config_options);
         exit_code = setup.run(slha_io);
         }
         break;
      case Gm2_cmd_line_options::THDM: {
         auto setup = make_thdm_setup(config_options);
         exit_code = setup.run(slha_io);
         }
         break;
      }
   } catch (const gm2calc::Error& error) {
      print_error(error, slha_io, config_options);
      exit_code = EXIT_FAILURE;
   }

   return exit_code;
}
