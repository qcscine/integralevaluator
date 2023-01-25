/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombConstructor.h>
#include <LibintIntegrals/TwoBodyIntegrals/SymmetryHelper.h>
#include <Utils/DataStructures/DensityMatrix.h>

namespace Scine {
namespace Integrals {
namespace TwoBody {

TwoTypeCoulombConstructor::TwoTypeCoulombConstructor(const Utils::DensityMatrix& densityMatrix_type1,
                                                     const Utils::DensityMatrix& densityMatrix_type2)
  : densityMatrix_type1_(densityMatrix_type1),
    densityMatrix_type2_(densityMatrix_type2),
    dim1_(densityMatrix_type1_.restrictedMatrix().rows()),
    dim2_(densityMatrix_type2_.restrictedMatrix().rows()) {
  coulomb_type1_.resizeLike(densityMatrix_type1_.restrictedMatrix());
  coulomb_type1_.setZero();
  coulomb_type2_.resizeLike(densityMatrix_type2_.restrictedMatrix());
  coulomb_type2_.setZero();
}

void TwoTypeCoulombConstructor::evaluateBasisFunctionQuartet(double integralValue, int basisFunction1,
                                                             int basisFunction2, int basisFunction3, int basisFunction4) {
  // Conveniently, the unrestricted and restricted cases can be handled together here.
  // The reason is that the restricted part of the unrestricted density matrix is the sum of the alpha and beta
  // density matrix.

  const auto* dm1 = densityMatrix_type1_.restrictedMatrix().data();
  const auto* dm2 = densityMatrix_type2_.restrictedMatrix().data();
  auto* J1 = coulomb_type1_.data();
  auto* J2 = coulomb_type2_.data();
  // Generates all non-equivalent permutations
  auto index = IndexType<4>({basisFunction1, basisFunction2, basisFunction3, basisFunction4});
  auto symmetricIndices = getSymmetricIndices<IntegralSymmetry::fourfold>(std::move(index));
  // Coulomb term
  for (const auto& idx : symmetricIndices) {
    J1[idx[0] * dim1_ + idx[1]] += dm2[idx[2] * dim2_ + idx[3]] * integralValue;
    J2[idx[2] * dim2_ + idx[3]] += dm1[idx[0] * dim1_ + idx[1]] * integralValue;
  }
}

void TwoTypeCoulombConstructor::finalizeEvaluation() {
  auto& L1 = coulomb_type1_;
  L1 = 0.5 * (L1 + L1.transpose()).eval();

  auto& L2 = coulomb_type2_;
  L2 = 0.5 * (L2 + L2.transpose()).eval();
}

const Eigen::MatrixXd& TwoTypeCoulombConstructor::getCoulombMatrixType1() const {
  return coulomb_type1_;
}

const Eigen::MatrixXd& TwoTypeCoulombConstructor::getCoulombMatrixType2() const {
  return coulomb_type2_;
}

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine
