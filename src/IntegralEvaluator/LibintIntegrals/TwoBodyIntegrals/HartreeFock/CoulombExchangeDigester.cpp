/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeDigester.h>

#define UNUSED(expr) \
  do {               \
    (void)(expr);    \
  } while (0)

namespace Scine {
namespace Integrals {
namespace TwoBody {

void CoulombExchangeDigester::digestImpl(double integralValue, int i, int j, int k, int l, int index, double degeneracy) {
  UNUSED(index);

  // IndexType<4> indexSet = {i , j, k, l};
  // if (!std::equal(indexSet.begin(), indexSet.end(), getMappedIndex<IntegralSymmetry::eightfold>(indexSet).begin())) {
  //  return;
  //}

  const auto threadNr = omp_get_thread_num();

  constructor_[threadNr].evaluateBasisFunctionQuartet(integralValue * this->scaling_ * degeneracy, i, j, k, l);
}

const std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix>& CoulombExchangeDigester::getResultImpl() const {
  return coulomb_exchange_;
}

double CoulombExchangeDigester::computeDegeneracyImpl(int shell1, int shell2, int shell3, int shell4) {
  auto shell12_deg = (shell1 == shell2) ? 1 : 2;
  auto shell34_deg = (shell3 == shell4) ? 1 : 2;
  auto shell12_34_deg = (shell1 == shell3) ? (shell2 == shell4 ? 1 : 2) : 2;
  return shell12_deg * shell34_deg * shell12_34_deg;
}

void CoulombExchangeDigester::initializeImpl(int numberThreads) {
  constructor_ = std::vector<CoulombExchangeConstructor>(numberThreads, CoulombExchangeConstructor(densityMatrix_));
}

void CoulombExchangeDigester::finalizeImpl() {
  for (auto& elem : constructor_) {
    elem.finalizeEvaluation();
  }

  if (densityMatrix_.restricted()) {
    for (auto& elem : constructor_) {
      coulomb_exchange_.first.restrictedMatrix() += elem.getCoulombMatrix().restrictedMatrix();
      coulomb_exchange_.second.restrictedMatrix() += elem.getExchangeMatrix().restrictedMatrix();
    }
  }
  else {
    for (auto& elem : constructor_) {
      coulomb_exchange_.first.alphaMatrix() += elem.getCoulombMatrix().alphaMatrix();
      coulomb_exchange_.second.alphaMatrix() += elem.getExchangeMatrix().alphaMatrix();
      if (densityMatrix_.numberElectronsInBetaMatrix() > 0) {
        coulomb_exchange_.first.betaMatrix() += elem.getCoulombMatrix().betaMatrix();
        coulomb_exchange_.second.betaMatrix() += elem.getExchangeMatrix().betaMatrix();
      }
    }
  }
}

CoulombExchangeDigester::CoulombExchangeDigester(const Utils::Integrals::BasisSet& scineBasis1,
                                                 const Utils::Integrals::BasisSet& scineBasis2,
                                                 const Utils::Integrals::IntegralSpecifier& specifier,
                                                 const Utils::DensityMatrix& densityMatrix)
  : Digester<CoulombExchangeDigester>(scineBasis1, scineBasis2, specifier), densityMatrix_(densityMatrix) {
  if (specifier.typeVector.size() == 2) {
    this->scaling_ = specifier.typeVector[0].charge * specifier.typeVector[1].charge;
  }
  // Derivatives not implemented, yet
  if (this->specifier_.derivOrder != 0) {
    throw std::runtime_error("Derivative of the Fock matrix not available, yet!");
  }

  // TODO see if this works properly:
  if (densityMatrix_.restricted()) {
    coulomb_exchange_.first.restrictedMatrix().resizeLike(densityMatrix_.restrictedMatrix());
    coulomb_exchange_.first.restrictedMatrix().setZero();
    coulomb_exchange_.second.restrictedMatrix().resizeLike(densityMatrix_.restrictedMatrix());
    coulomb_exchange_.second.restrictedMatrix().setZero();
  }
  else {
    coulomb_exchange_.first.alphaMatrix().resizeLike(densityMatrix_.alphaMatrix());
    coulomb_exchange_.first.alphaMatrix().setZero();
    coulomb_exchange_.second.alphaMatrix().resizeLike(densityMatrix_.alphaMatrix());
    coulomb_exchange_.second.alphaMatrix().setZero();
    // TODO see if this works:
    if (densityMatrix_.numberElectronsInBetaMatrix() > 0) {
      coulomb_exchange_.first.betaMatrix().resizeLike(densityMatrix_.betaMatrix());
      coulomb_exchange_.first.betaMatrix().setZero();
      coulomb_exchange_.second.betaMatrix().resizeLike(densityMatrix_.betaMatrix());
      coulomb_exchange_.second.betaMatrix().setZero();
    }
  }
}

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine
