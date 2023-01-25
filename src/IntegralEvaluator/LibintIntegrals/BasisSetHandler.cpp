/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* internal */
#include <LibintIntegrals/BasisSetHandler.h>
#include <Utils/DataStructures/IntegralSpecifier.h>
#include <Utils/Geometry/ElementInfo.h>
/* external */
#include <libint2/util/small_vector.h>

using namespace Scine;
using namespace Integrals;

auto BasisSetHandler::libintToScine(const libint2::BasisSet& libintBasis, const Utils::AtomCollection& atoms,
                                    const bool& pure_solid) -> Utils::Integrals::BasisSet {
  Utils::Integrals::BasisSet scineBasis(atoms);

  for (const auto& libintShell : libintBasis) {
    Utils::Displacement shift;
    std::size_t max_l;

    // alpha is same for all contractions
    std::vector<double> v_alpha;
    std::copy(libintShell.alpha.begin(), libintShell.alpha.end(), std::back_inserter(v_alpha));
    // shift as well
    for (int i = 0; i < 3; ++i) {
      shift(i) = libintShell.O.at(i);
    }
    // Loop over contractions
    // Remember: a libint shell can have multiple contracted orbitals, but a scine shell only one.
    for (auto const& contr : libintShell.contr) {
      assert(contr.pure == pure_solid);
      // maximum angular momentum
      max_l = contr.l;
      // coefficients are different in contractions
      std::vector<double> v_coeffs;
      std::copy(contr.coeff.begin(), contr.coeff.end(), std::back_inserter(v_coeffs));
      Utils::Integrals::Shell new_shell(v_alpha, v_coeffs, shift, max_l, pure_solid);
      scineBasis.emplace_back(new_shell);
    }
  }

  scineBasis.setPureSpherical(pure_solid);

  return scineBasis;
}

auto BasisSetHandler::scineToLibint(const Utils::AtomCollection& scineAtoms) -> std::vector<libint2::Atom> {
  std::vector<libint2::Atom> libintAtoms;
  libintAtoms.reserve(scineAtoms.size());

  for (const auto& atom : scineAtoms) {
    libintAtoms.push_back({static_cast<int>(Utils::ElementInfo::Z(atom.getElementType())), atom.getPosition()(0),
                           atom.getPosition()(1), atom.getPosition()(2)});
  }

  return libintAtoms;
}

auto BasisSetHandler::scineToLibint(const Utils::Integrals::Operator& scineOp) -> libint2::Operator {
  switch (scineOp) {
    case Utils::Integrals::Operator::Overlap:
      return libint2::Operator::overlap;
    case Utils::Integrals::Operator::Kinetic:
      return libint2::Operator::kinetic;
    case Utils::Integrals::Operator::PointCharges:
      return libint2::Operator::nuclear;
    case Utils::Integrals::Operator::Coulomb:
      return libint2::Operator::coulomb;
    case Utils::Integrals::Operator::CoulombCOM:
      throw std::runtime_error("Conversion from CoulombCOM to libint2 operator not possible!");
    case Utils::Integrals::Operator::Dipole:
      return libint2::Operator::emultipole1;
    case Utils::Integrals::Operator::KineticCOM:
      return libint2::Operator::kinetic;
    default:
      throw std::runtime_error("Operator not recognized!");
  }
}

auto BasisSetHandler::generateShellPairs(Utils::Integrals::BasisSet& basis, bool performOverlapPrescreening,
                                         double threshold, bool calculateCauchySchwarzFactor) -> void {
  Libint::getInstance();
  auto shellPairs = std::make_shared<Utils::Integrals::ShellPairs>();
  shellPairs->resize(basis.size());
#pragma omp parallel
  {
    auto localEngine = Libint::getEngine(basis, libint2::Operator::overlap);
    auto const& buffer = localEngine.results();
#pragma omp barrier
#pragma omp for schedule(dynamic)
    for (size_t s1 = 0; s1 < basis.size(); ++s1) {
      const auto& shell1 = basis[s1];
      Utils::Integrals::ShellPair shellPair(shell1, calculateCauchySchwarzFactor);

      auto const shell1Size = basis[s1].size();
      for (size_t s2 = 0; s2 <= s1; ++s2) {
        auto const shell2Size = basis[s2].size();
        const auto& shell2 = basis[s2];

        if (performOverlapPrescreening) {
          bool onSameCenter = shell1.getShift() == shell2.getShift();
          bool toBeIncluded = onSameCenter;

          if (!onSameCenter) {
            localEngine.compute(scineToLibint(shell1), scineToLibint(shell2));
            double overlapNorm =
                (Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(buffer[0], shell1Size, shell2Size)).norm();
            toBeIncluded = overlapNorm > threshold;
          }
          if (toBeIncluded) {
            addPair(shellPair, s2, shell2);
          }
        }
        else {
          addPair(shellPair, s2, shell2);
        }
      }
      std::sort(shellPair.begin(), shellPair.end(),
                [&](const Utils::Integrals::ShellPairData& first, const Utils::Integrals::ShellPairData& second) {
                  return first.secondShellIndex < second.secondShellIndex;
                });
#pragma omp critical(shellPairAddition)
      { shellPairs->at(s1) = std::move(shellPair); }
    }
  }
  shellPairs->setCauchySchwarzFactor(calculateCauchySchwarzFactor);

  basis.setShellPairs(shellPairs);
}
void BasisSetHandler::addPair(Utils::Integrals::ShellPair& ShellPair, size_t secondShell, const Utils::Integrals::Shell& shell2,
                              double ln_prec, const bool calculateCauchySchwarzFactor) {
  Utils::Integrals::ShellPairData newPair;
  newPair.secondShellIndex = secondShell;
  const auto libintShell1 = scineToLibint(ShellPair.getShell());
  const auto libintShell2 = scineToLibint(shell2);
  newPair.precomputedShellPair = std::make_unique<Utils::Integrals::ShellPairType>(ShellPair.getShell(), shell2, ln_prec);

  if (calculateCauchySchwarzFactor) {
    Libint::getInstance();

    auto const shell1Size = libintShell1.size();
    auto const shell2Size = libintShell2.size();
    libint2::Engine localEngine;
#pragma omp critical(initializeLocalEngine)
    {
      localEngine = Libint::getEngine(libint2::Operator::coulomb, libintShell1, libintShell2, libintShell1, libintShell2);
    }

    // Important for Cauchy-Schwarz that no native screening is performed!!
    localEngine.set_precision(0);
    auto const& buffer = localEngine.results();
    localEngine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(libintShell1, libintShell2,
                                                                                libintShell1, libintShell2);
    newPair.cauchySchwarzFactor =
        std::sqrt(Eigen::Map<const Eigen::MatrixXd>(buffer[0], shell1Size, shell2Size).cwiseAbs().maxCoeff());
  }
#pragma omp critical(addDataToShell)
  { ShellPair.emplace_back(std::move(newPair)); }
}
