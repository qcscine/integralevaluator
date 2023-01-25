/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_BASISSETHANDLER_H
#define INTEGRALEVALUATOR_BASISSETHANDLER_H

#include <LibintIntegrals/Libint.h>
#include <Utils/DataStructures/BasisSet.h>
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace Utils {
namespace Integrals {
class ShellPair;
class ShellPairs;
enum class Operator;
} // namespace Integrals
} // namespace Utils
} // namespace Scine

namespace Scine {
namespace Integrals {

class BasisSetHandler {
 public:
  BasisSetHandler() = default;
  ~BasisSetHandler() = default;

  /**
   * @brief converts a `libint2::BasisSet` object to a `Utils::Integrals::BasisSet` object.
   * @param libintBasis
   * @param atoms
   * @param pure_solid
   * @return scineBasis
   */
  static auto libintToScine(const libint2::BasisSet& libintBasis, const Utils::AtomCollection& atoms, const bool& pure_solid)
      -> Utils::Integrals::BasisSet;

  /**
   * @brief converts Scine `AtomCollection` to `vector<libint2::Atom>`
   * @param scineAtoms
   * @return vector of libint2 atoms.
   */
  static auto scineToLibint(const Utils::AtomCollection& scineAtoms) -> std::vector<libint2::Atom>;

  /**
   * @brief Converts `scine` to `libint` operator enum
   * @param scine operator
   * @return libint2::Operator
   */
  static auto scineToLibint(const Utils::Integrals::Operator& scineOp) -> libint2::Operator;
  /**
   * @brief Converts a Scine `Shell` to a Libint `Shell`.
   * @param scine Shell.
   * @return libint Shell.
   */
  inline static auto scineToLibint(const Utils::Integrals::Shell& scineShell) -> libint2::Shell {
    libint2::Shell libintShell;

    libintShell.alpha.resize(scineShell.getVecAlpha().size());
    for (std::size_t it = 0; it < scineShell.getVecAlpha().size(); ++it) {
      libintShell.alpha.at(it) = scineShell.getVecAlpha().at(it);
    }

    libintShell.max_ln_coeff.resize(scineShell.getMaxLnCoeffs().size());
    for (std::size_t it = 0; it < scineShell.getMaxLnCoeffs().size(); ++it)
      libintShell.max_ln_coeff.at(it) = scineShell.getMaxLnCoeffs().at(it);

    libint2::svector<double> libintCoeffs(scineShell.getVecCoeffs().size());
    for (std::size_t it = 0; it < scineShell.getVecCoeffs().size(); ++it)
      libintCoeffs.at(it) = scineShell.getVecCoeffs().at(it);

    libintShell.contr = {{static_cast<int>(scineShell.l()), scineShell.isPureSolid(), libintCoeffs}};
    libintShell.O = {{scineShell.getShift().data()[0], scineShell.getShift().data()[1], scineShell.getShift().data()[2]}};

    return libintShell;
  }
  /**
   * @brief Converts a Scine `ShellPairType` object to a Libint `ShellPair`.
   * note: is made inline because it is required in a performance intensive part.
   * @param scine ShellPairType
   * @return libint Shell Pair
   */
  inline static auto scineToLibint(const Utils::Integrals::ShellPairType& scineShellPairType) -> libint2::ShellPair {
    libint2::ShellPair libintShellPair;

    libintShellPair.primpairs.reserve(scineShellPairType.primpairs.size());
    libintShellPair.screening_method_ = libint2::ScreeningMethod::Original;

    for (auto i = 0UL; i < scineShellPairType.primpairs.size(); ++i) {
      libint2::ShellPair::PrimPairData data;
      std::copy(std::begin(scineShellPairType.primpairs[i].P), std::end(scineShellPairType.primpairs[i].P),
                std::begin(data.P));
      data.K = scineShellPairType.primpairs[i].K;
      data.one_over_gamma = scineShellPairType.primpairs[i].one_over_gamma;
      data.ln_scr = scineShellPairType.primpairs[i].scr;
      data.p1 = scineShellPairType.primpairs[i].p1;
      data.p2 = scineShellPairType.primpairs[i].p2;
      libintShellPair.primpairs.push_back(data);
    }
    std::copy(std::begin(scineShellPairType.AB), std::end(scineShellPairType.AB), std::begin(libintShellPair.AB));
    return libintShellPair;
  }
  /**
   * @brief Adds a shell interaction to a given shell.
   * @param ln_prec This is taken as suggested by the libint hartree-fock++.cc example file.
   */
  static void addPair(Utils::Integrals::ShellPair& ShellPair, size_t secondShell, const Utils::Integrals::Shell& shell2,
                      double ln_prec = std::log(std::numeric_limits<double>::epsilon() / 1e10),
                      bool calculateCauchySchwarzFactor = true);

  static auto generateShellPairs(Utils::Integrals::BasisSet& scineBasis, bool performOverlapPrescreening = true,
                                 double threshold = 1e-12, bool calculateCauchySchwarzFactor = true) -> void;
};

} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_BASISSETHANDLER_H
