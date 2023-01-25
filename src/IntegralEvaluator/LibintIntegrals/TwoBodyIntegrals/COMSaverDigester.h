/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_COMSAVERDIGESTER_H
#define INTEGRALEVALUATOR_COMSAVERDIGESTER_H

#include <LibintIntegrals/LibintIntegrals.h>
#include <LibintIntegrals/TwoBodyIntegrals/Digester.h>
#include <LibintIntegrals/TwoBodyIntegrals/SymmetryHelper.h>
#include <Eigen/Dense>

namespace Scine {
namespace Integrals {
namespace TwoBody {

template<IntegralSymmetry symmetry>
class COMSaverDigester : public Digester<COMSaverDigester<symmetry>> {
 public:
  COMSaverDigester(const Utils::Integrals::BasisSet& scineBasis1, const Utils::Integrals::BasisSet& scineBasis2,
                   const Utils::Integrals::IntegralSpecifier& specifier);

  void digestImpl(double integralValue, int basisFunction1, int basisFunction2, int basisFunction3, int basisFunction4,
                  int index, double degeneracy);

  double computeDegeneracyImpl(int shell1, int shell2, int shell3, int shell4);

  const IntegralEvaluatorMap& getResultImpl() const;
  void initializeImpl(int numberThreads);
  void finalizeImpl();

 private:
  IntegralEvaluatorMap result_;

  std::vector<double*> resultPtr_;

  std::vector<IntegralEvaluatorMap> vectorMomentumIntegrals_;
  double totalMass_;
};

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_COMSAVERDIGESTER_H
