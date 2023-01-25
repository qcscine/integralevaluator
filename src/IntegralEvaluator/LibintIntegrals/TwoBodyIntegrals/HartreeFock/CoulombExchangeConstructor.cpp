/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeConstructor.h>
#include <Utils/DataStructures/DensityMatrix.h>

namespace Scine {
namespace Integrals {
namespace TwoBody {

CoulombExchangeConstructor::CoulombExchangeConstructor(const Utils::DensityMatrix& densityMatrix)
  : densityMatrix_(densityMatrix), dim_(densityMatrix_.restrictedMatrix().rows()) {
  if (densityMatrix_.restricted()) {
    exchange_.restrictedMatrix().resizeLike(densityMatrix_.restrictedMatrix());
    exchange_.restrictedMatrix().setZero();
    coulomb_.restrictedMatrix().resizeLike(densityMatrix_.restrictedMatrix());
    coulomb_.restrictedMatrix().setZero();
  }
  else {
    exchange_.alphaMatrix().resizeLike(densityMatrix_.alphaMatrix());
    exchange_.alphaMatrix().setZero();
    coulomb_.alphaMatrix().resizeLike(densityMatrix_.alphaMatrix());
    coulomb_.alphaMatrix().setZero();
    // This relies on the density matrix being built properly.
    if (densityMatrix_.numberElectronsInBetaMatrix() > 0) {
      exchange_.betaMatrix().resizeLike(densityMatrix_.betaMatrix());
      exchange_.betaMatrix().setZero();
      coulomb_.betaMatrix().resizeLike(densityMatrix_.betaMatrix());
      coulomb_.betaMatrix().setZero();
    }
  }
}

void CoulombExchangeConstructor::evaluateBasisFunctionQuartet(double integralValue, int basisFunction1,
                                                              int basisFunction2, int basisFunction3, int basisFunction4) {
  auto integralValue_05 = integralValue * 0.5;
  auto integralValue_025 = integralValue * 0.25;
  auto b1dim = basisFunction1 * dim_;
  auto b2dim = basisFunction2 * dim_;
  auto b3dim = basisFunction3 * dim_;
  // auto b4dim=basisFunction4*dim_;

  if (densityMatrix_.restricted()) {
    auto* J = coulomb_.restrictedMatrix().data();
    auto* K = exchange_.restrictedMatrix().data();
    const auto* dm = densityMatrix_.restrictedMatrix().data();

    J[b1dim + basisFunction2] += dm[b3dim + basisFunction4] * integralValue_05;
    J[b3dim + basisFunction4] += dm[b1dim + basisFunction2] * integralValue_05;
    K[b1dim + basisFunction3] += dm[b2dim + basisFunction4] * integralValue_025;
    K[b2dim + basisFunction4] += dm[b1dim + basisFunction3] * integralValue_025;
    K[b1dim + basisFunction4] += dm[b2dim + basisFunction3] * integralValue_025;
    K[b2dim + basisFunction3] += dm[b1dim + basisFunction4] * integralValue_025;
  }
  else {
    auto* J_alpha = coulomb_.alphaMatrix().data();
    auto* K_alpha = exchange_.alphaMatrix().data();
    const auto* dm_alpha = densityMatrix_.alphaMatrix().data();

    J_alpha[b1dim + basisFunction2] += dm_alpha[b3dim + basisFunction4] * integralValue_05;
    J_alpha[b3dim + basisFunction4] += dm_alpha[b1dim + basisFunction2] * integralValue_05;
    K_alpha[b1dim + basisFunction3] += dm_alpha[b2dim + basisFunction4] * integralValue_025;
    K_alpha[b2dim + basisFunction4] += dm_alpha[b1dim + basisFunction3] * integralValue_025;
    K_alpha[b1dim + basisFunction4] += dm_alpha[b2dim + basisFunction3] * integralValue_025;
    K_alpha[b2dim + basisFunction3] += dm_alpha[b1dim + basisFunction4] * integralValue_025;

    if (densityMatrix_.numberElectronsInBetaMatrix() > 0) {
      auto* J_beta = coulomb_.betaMatrix().data();
      auto* K_beta = exchange_.betaMatrix().data();
      const auto* dm_beta = densityMatrix_.betaMatrix().data();

      J_beta[b1dim + basisFunction2] += dm_beta[b3dim + basisFunction4] * integralValue_05;
      J_beta[b3dim + basisFunction4] += dm_beta[b1dim + basisFunction2] * integralValue_05;
      K_beta[b1dim + basisFunction3] += dm_beta[b2dim + basisFunction4] * integralValue_025;
      K_beta[b2dim + basisFunction4] += dm_beta[b1dim + basisFunction3] * integralValue_025;
      K_beta[b1dim + basisFunction4] += dm_beta[b2dim + basisFunction3] * integralValue_025;
      K_beta[b2dim + basisFunction3] += dm_beta[b1dim + basisFunction4] * integralValue_025;
    }
  }
}

void CoulombExchangeConstructor::finalizeEvaluation() {
  if (!densityMatrix_.restricted()) {
    auto& J_alpha = coulomb_.alphaMatrix();
    auto& K_alpha = exchange_.alphaMatrix();
    J_alpha = 0.5 * (J_alpha + J_alpha.transpose()).eval();
    K_alpha = 0.5 * (K_alpha + K_alpha.transpose()).eval();

    if (densityMatrix_.numberElectronsInBetaMatrix() > 0) {
      auto& J_beta = coulomb_.betaMatrix();
      auto& K_beta = exchange_.betaMatrix();
      J_beta = 0.5 * (J_beta + J_beta.transpose()).eval();
      K_beta = 0.5 * (K_beta + K_beta.transpose()).eval();
    }
  }
  else {
    auto& J_restricted = coulomb_.restrictedMatrix();
    auto& K_restricted = exchange_.restrictedMatrix();
    J_restricted = 0.5 * (J_restricted + J_restricted.transpose()).eval();
    K_restricted = 0.5 * (K_restricted + K_restricted.transpose()).eval();
  }
}

const Utils::SpinAdaptedMatrix& CoulombExchangeConstructor::getCoulombMatrix() const {
  return coulomb_;
}

const Utils::SpinAdaptedMatrix& CoulombExchangeConstructor::getExchangeMatrix() const {
  return exchange_;
}

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine
