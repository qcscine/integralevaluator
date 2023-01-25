/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <LibintIntegrals/Libint.h>
#include <LibintIntegrals/TwoBodyIntegrals/CauchySchwarzDensityPrescreener.h>
#include <Utils/DataStructures/BasisSet.h>
#include <Utils/DataStructures/DensityMatrix.h>

using namespace Scine;
using namespace Integrals;
using namespace TwoBody;

CauchySchwarzDensityPrescreener::CauchySchwarzDensityPrescreener(const Utils::Integrals::BasisSet& basisSet,
                                                                 const Utils::DensityMatrix& densityMatrix,
                                                                 double prescreenThreshold)
  : TwoBodiesIntegralsPrescreener<CauchySchwarzDensityPrescreener>(basisSet, prescreenThreshold),
    densityMatrix_(densityMatrix) {
  calculateShellBlockDensityMatrix();
}

bool CauchySchwarzDensityPrescreener::isSignificantImpl(int shell1, int shell2, int shell3, int shell4,
                                                        double cauchySchwarzFactor) const {
  // if no Cauchy-Schwarz factors were calculated, do not perform screening.
  if (!scineBasis_.getShellPairs()->hasCauchySchwarzFactor())
    return true;

  // Very easy to check, do screening with biggest element of density matrix.
  if (std::abs(densityMaximum_ * cauchySchwarzFactor) < prescreeningThreshold_) {
    return false;
  }
  double maxDensityShellBlock =
      std::max(densityMatrixShellBlockMaxima_(shell1, shell2), densityMatrixShellBlockMaxima_(shell1, shell3));
  maxDensityShellBlock = std::max(maxDensityShellBlock, densityMatrixShellBlockMaxima_(shell1, shell4));
  maxDensityShellBlock = std::max(maxDensityShellBlock, densityMatrixShellBlockMaxima_(shell2, shell3));
  maxDensityShellBlock = std::max(maxDensityShellBlock, densityMatrixShellBlockMaxima_(shell2, shell4));
  maxDensityShellBlock = std::max(maxDensityShellBlock, densityMatrixShellBlockMaxima_(shell3, shell4));
  return (std::abs(maxDensityShellBlock * cauchySchwarzFactor) > prescreeningThreshold_);
}

void CauchySchwarzDensityPrescreener::calculateShellBlockDensityMatrix() {
  auto const& s2bf = scineBasis_.shell2bf();
  densityMatrixShellBlockMaxima_.resize(scineBasis_.size(), scineBasis_.size());

  for (auto shell1 = 0UL; shell1 < scineBasis_.size(); ++shell1) {
    for (auto shell2 = 0UL; shell2 <= shell1; ++shell2) {
      const Eigen::MatrixXd shellBlock = densityMatrix_.restrictedMatrix().block(
          s2bf[shell1], s2bf[shell2], scineBasis_[shell1].size(), scineBasis_[shell2].size());
      densityMatrixShellBlockMaxima_(shell1, shell2) = shellBlock.cwiseAbs().maxCoeff();
    }
  }
  densityMaximum_ = densityMatrixShellBlockMaxima_.maxCoeff();
  densityMatrixShellBlockMaxima_ = densityMatrixShellBlockMaxima_.selfadjointView<Eigen::Lower>();
}
