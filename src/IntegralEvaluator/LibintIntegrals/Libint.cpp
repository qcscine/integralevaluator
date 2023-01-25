/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <LibintIntegrals/Libint.h>

using namespace Scine;
using namespace Integrals;

int Libint::nThreads_ = 1;

Libint::Libint() {
  libint2::initialize();
  nThreads_ = omp_get_max_threads();
}

Libint::~Libint() {
  libint2::finalize();
}

auto Libint::getEngine(const Scine::Utils::Integrals::BasisSet& basis1, const Scine::Utils::Integrals::BasisSet& basis2,
                       const libint2::Operator& op, const int& derivOrder) -> libint2::Engine {
  getInstance();
  libint2::Engine engine;
  return {op, std::max(basis1.max_nprim(), basis2.max_nprim()),
          static_cast<int>(std::max(basis1.max_l(), basis2.max_l())), derivOrder};
}

auto Libint::getEngine(const libint2::Operator& op, const libint2::Shell& shell1, const libint2::Shell& shell2,
                       const libint2::Shell& shell3, const libint2::Shell& shell4) -> libint2::Engine {
  getInstance();
  int angularMomentum = 0;
  for (const auto& primitive : shell1.contr) {
    angularMomentum = std::max(primitive.l, angularMomentum);
  }
  for (const auto& primitive : shell2.contr) {
    angularMomentum = std::max(primitive.l, angularMomentum);
  }
  for (const auto& primitive : shell3.contr) {
    angularMomentum = std::max(primitive.l, angularMomentum);
  }
  for (const auto& primitive : shell4.contr) {
    angularMomentum = std::max(primitive.l, angularMomentum);
  }
  auto nPrimitives = shell1.nprim();
  nPrimitives = std::max(shell2.nprim(), nPrimitives);
  nPrimitives = std::max(shell3.nprim(), nPrimitives);
  nPrimitives = std::max(shell4.nprim(), nPrimitives);

  return {op, nPrimitives, angularMomentum};
}

auto Libint::getEngine(const Scine::Utils::Integrals::BasisSet& basis, const libint2::Operator& op, const int& derivOrder)
    -> libint2::Engine {
  getInstance();
  return {op, basis.max_nprim(), static_cast<int>(basis.max_l()), derivOrder};
}
