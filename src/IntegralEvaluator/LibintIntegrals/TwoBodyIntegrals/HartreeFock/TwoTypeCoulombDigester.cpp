/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombDigester.h>

#define UNUSED(expr) \
  do {               \
    (void)(expr);    \
  } while (0)

namespace Scine {
namespace Integrals {
namespace TwoBody {

void TwoTypeCoulombDigester::digestImpl(double integralValue, int i, int j, int k, int l, int index, double degeneracy) {
  UNUSED(degeneracy);
  UNUSED(index);

  integralValue *= this->scaling_;

  IndexType<4> indexSet = {i, j, k, l};
  if (!std::equal(indexSet.begin(), indexSet.end(), getMappedIndex<IntegralSymmetry::fourfold>(indexSet).begin())) {
    return;
  }

  const auto threadNr = omp_get_thread_num();

  constructor_[threadNr].evaluateBasisFunctionQuartet(integralValue, i, j, k, l);
}

const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& TwoTypeCoulombDigester::getResultImpl() const {
  return coulomb_type1_type2_;
}

double TwoTypeCoulombDigester::computeDegeneracyImpl(int shell1, int shell2, int shell3, int shell4) {
  UNUSED(shell1);
  UNUSED(shell2);
  UNUSED(shell3);
  UNUSED(shell4);
  return 1;
}

void TwoTypeCoulombDigester::initializeImpl(int numberThreads) {
  constructor_ = std::vector<TwoTypeCoulombConstructor>(
      numberThreads, TwoTypeCoulombConstructor(densityMatrix_type1_, densityMatrix_type2_));
}

void TwoTypeCoulombDigester::finalizeImpl() {
  for (auto& elem : constructor_) {
    elem.finalizeEvaluation();
  }

  for (auto& elem : constructor_) {
    coulomb_type1_type2_.first += elem.getCoulombMatrixType1();
  }
  for (auto& elem : constructor_) {
    coulomb_type1_type2_.second += elem.getCoulombMatrixType2();
  }
}

TwoTypeCoulombDigester::TwoTypeCoulombDigester(const Utils::Integrals::BasisSet& scineBasis1,
                                               const Utils::Integrals::BasisSet& scineBasis2,
                                               const Utils::Integrals::IntegralSpecifier& specifier,
                                               const Utils::DensityMatrix& densityMatrix1,
                                               const Utils::DensityMatrix& densityMatrix2)
  : Digester<TwoTypeCoulombDigester>(scineBasis1, scineBasis2, specifier),
    densityMatrix_type1_(densityMatrix1),
    densityMatrix_type2_(densityMatrix2) {
  if (specifier.typeVector.size() == 2) {
    this->scaling_ = specifier.typeVector[0].charge * specifier.typeVector[1].charge;
  }

  // Derivatives not implemented, yet
  if (this->specifier_.derivOrder != 0) {
    throw std::runtime_error("Derivative of the Fock matrix not available, yet!");
  }

  coulomb_type1_type2_.first.resizeLike(densityMatrix_type1_.restrictedMatrix());
  coulomb_type1_type2_.first.setZero();
  coulomb_type1_type2_.second.resizeLike(densityMatrix_type2_.restrictedMatrix());
  coulomb_type1_type2_.second.setZero();
}

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine
