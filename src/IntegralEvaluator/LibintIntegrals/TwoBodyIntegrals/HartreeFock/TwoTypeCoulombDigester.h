/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INTEGRALEVALUATOR_TWOTYPECOULOMBDIGESTER_H
#define INTEGRALEVALUATOR_TWOTYPECOULOMBDIGESTER_H

#include <LibintIntegrals/LibintIntegrals.h>
#include <LibintIntegrals/TwoBodyIntegrals/Digester.h>
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombConstructor.h>
#include <LibintIntegrals/TwoBodyIntegrals/SymmetryHelper.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Eigen/Dense>

namespace Scine {
namespace Integrals {
namespace TwoBody {

class TwoTypeCoulombDigester : public Digester<TwoTypeCoulombDigester> {
 public:
  TwoTypeCoulombDigester(const Utils::Integrals::BasisSet& scineBasis1, const Utils::Integrals::BasisSet& scineBasis2,
                         const Utils::Integrals::IntegralSpecifier& specifier, const Utils::DensityMatrix& D1,
                         const Utils::DensityMatrix& D2);

  void digestImpl(double integralValue, int basisFunction1, int basisFunction2, int basisFunction3, int basisFunction4,
                  int index, double degeneracy);

  double computeDegeneracyImpl(int shell1, int shell2, int shell3, int shell4);

  const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& getResultImpl() const;
  void initializeImpl(int numberThreads);
  void finalizeImpl();

 private:
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> coulomb_type1_type2_;
  const Utils::DensityMatrix& densityMatrix_type1_;
  const Utils::DensityMatrix& densityMatrix_type2_;
  /* One constructor per thread */
  std::vector<TwoTypeCoulombConstructor> constructor_;
};

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_TWOTYPECOULOMBDIGESTER_H
