/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_CAUCHYSCHWARZDENSITYPRESCREENER_H
#define INTEGRALEVALUATOR_CAUCHYSCHWARZDENSITYPRESCREENER_H

#include <LibintIntegrals/TwoBodyIntegrals/Prescreener.h>
#include <Eigen/Core>

namespace Scine {
namespace Integrals {
namespace TwoBody {

class CauchySchwarzDensityPrescreener : public TwoBodiesIntegralsPrescreener<CauchySchwarzDensityPrescreener> {
 public:
  CauchySchwarzDensityPrescreener(const Utils::Integrals::BasisSet& basisSet, const Utils::DensityMatrix& densityMatrix,
                                  double prescreenThreshold = 1e-12);

  bool isSignificantImpl(int shell1, int shell2, int shell3, int shell4, double preescreenThreshold) const;

 private:
  void calculateShellBlockDensityMatrix();
  const Utils::DensityMatrix& densityMatrix_;
  Eigen::MatrixXd densityMatrixShellBlockMaxima_;
  double densityMaximum_{};
};

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_CAUCHYSCHWARZDENSITYPRESCREENER_H
