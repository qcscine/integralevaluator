/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INTEGRALEVALUATOR_COULOMBEXCHANGEDIGESTER_H
#define INTEGRALEVALUATOR_COULOMBEXCHANGEDIGESTER_H

#include <LibintIntegrals/LibintIntegrals.h>
#include <LibintIntegrals/TwoBodyIntegrals/Digester.h>
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeConstructor.h>
#include <LibintIntegrals/TwoBodyIntegrals/SymmetryHelper.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Eigen/Dense>

namespace Scine {
namespace Integrals {
namespace TwoBody {

class CoulombExchangeDigester : public Digester<CoulombExchangeDigester> {
 public:
  CoulombExchangeDigester(const Utils::Integrals::BasisSet& scineBasis1, const Utils::Integrals::BasisSet& scineBasis2,
                          const Utils::Integrals::IntegralSpecifier& specifier, const Utils::DensityMatrix& D);

  void digestImpl(double integralValue, int basisFunction1, int basisFunction2, int basisFunction3, int basisFunction4,
                  int index, double degeneracy);

  double computeDegeneracyImpl(int shell1, int shell2, int shell3, int shell4);

  const std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix>& getResultImpl() const;
  void initializeImpl(int numberThreads);
  void finalizeImpl();

 private:
  std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix> coulomb_exchange_;
  const Utils::DensityMatrix& densityMatrix_;
  /* One constructor per thread */
  std::vector<CoulombExchangeConstructor> constructor_;
};

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_COULOMBEXCHANGEDIGESTER_H
