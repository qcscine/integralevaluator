/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <LibintIntegrals/OneBodyIntegrals.h>
#include <LibintIntegrals/TwoBodyIntegrals/COMSaverDigester.h>

#define UNUSED(expr) \
  do {               \
    (void)(expr);    \
  } while (0)

namespace Scine {
namespace Integrals {
namespace TwoBody {

template<>
void COMSaverDigester<IntegralSymmetry::fourfold>::digestImpl(double integralValue, int i, int j, int k, int l,
                                                              int index, double degeneracy) {
  UNUSED(degeneracy);

  auto& ptr_data = resultPtr_[index];

  double COM = 0;
  for (auto const& derivKey : {Utils::Integrals::DerivKey::x, Utils::Integrals::DerivKey::y, Utils::Integrals::DerivKey::z}) {
    COM -= (1 / totalMass_) * vectorMomentumIntegrals_[0][{Utils::Integrals::Component::none, derivKey, 0}](i, j) *
           vectorMomentumIntegrals_[0][{Utils::Integrals::Component::none, derivKey, 0}](k, l);
  }
  integralValue *= this->scaling_;
  // two-fold
  // Subtract COM:
  // ij,kl = ji,lk
  ptr_data[i * dim1_ + j + dim1sq_ * (k * dim2_ + l)] = integralValue - COM;
  ptr_data[j * dim1_ + i + dim1sq_ * (l * dim2_ + k)] = integralValue - COM;

  // Add COM:
  // ji,kl = ij,lk
  ptr_data[j * dim1_ + i + dim1sq_ * (k * dim2_ + l)] = integralValue + COM;
  ptr_data[i * dim1_ + j + dim1sq_ * (l * dim2_ + k)] = integralValue + COM;

  //
  // four-fold
  //
  // Subtract COM:
  // kl,ij =  = lk,ji
  ptr_data[k * dim1_ + l + dim1sq_ * (i * dim2_ + j)] = integralValue - COM;
  ptr_data[l * dim1_ + k + dim1sq_ * (j * dim2_ + i)] = integralValue - COM;
  // Add COM:
  // lk,ij  = kl,ji
  ptr_data[l * dim1_ + k + dim1sq_ * (i * dim2_ + j)] = integralValue + COM;
  ptr_data[k * dim1_ + l + dim1sq_ * (j * dim2_ + i)] = integralValue + COM;
}

template<>
void COMSaverDigester<IntegralSymmetry::twofold>::digestImpl(double integralValue, int i, int j, int k, int l,
                                                             int index, double degeneracy) {
  UNUSED(degeneracy);

  auto& ptr_data = resultPtr_[index];

  double COM = 0;

  for (auto const& derivKey : {Utils::Integrals::DerivKey::x, Utils::Integrals::DerivKey::y, Utils::Integrals::DerivKey::z}) {
    COM -= (1 / totalMass_) * vectorMomentumIntegrals_[0][{Utils::Integrals::Component::none, derivKey, 0}](i, j) *
           vectorMomentumIntegrals_[1][{Utils::Integrals::Component::none, derivKey, 0}](k, l);
  }
  integralValue *= this->scaling_;
  // two-fold
  // Subtract COM:
  // ij,kl = ji,lk
  ptr_data[i * dim1_ + j + dim1sq_ * (k * dim2_ + l)] = integralValue - COM;
  ptr_data[j * dim1_ + i + dim1sq_ * (l * dim2_ + k)] = integralValue - COM;

  // Add COM:
  // ji,kl = ij,lk
  ptr_data[j * dim1_ + i + dim1sq_ * (k * dim2_ + l)] = integralValue + COM;
  ptr_data[i * dim1_ + j + dim1sq_ * (l * dim2_ + k)] = integralValue + COM;
}

template<IntegralSymmetry symmetry>
const IntegralEvaluatorMap& COMSaverDigester<symmetry>::getResultImpl() const {
  return result_;
}

template<IntegralSymmetry symmetry>
double COMSaverDigester<symmetry>::computeDegeneracyImpl(int shell1, int shell2, int shell3, int shell4) {
  UNUSED(shell1);
  UNUSED(shell2);
  UNUSED(shell3);
  UNUSED(shell4);
  return 1;
}

template<IntegralSymmetry symmetry>
void COMSaverDigester<symmetry>::initializeImpl(int numberThreads) {
  UNUSED(numberThreads);
}

template<IntegralSymmetry symmetry>
void COMSaverDigester<symmetry>::finalizeImpl() {
}

template<IntegralSymmetry symmetry>
COMSaverDigester<symmetry>::COMSaverDigester(const Utils::Integrals::BasisSet& scineBasis1,
                                             const Utils::Integrals::BasisSet& scineBasis2,
                                             const Utils::Integrals::IntegralSpecifier& specifier)
  : Digester<COMSaverDigester<symmetry>>(scineBasis1, scineBasis2, specifier) {
  static_assert(symmetry == IntegralSymmetry::fourfold || symmetry == IntegralSymmetry::twofold,
                "Symmetry must be either 2- or 4-fold.");

  if (specifier.typeVector.size() == 2) {
    this->scaling_ = specifier.typeVector[0].charge * specifier.typeVector[1].charge;
  }
  for (int center = 0; center < this->numberOfCenters_; ++center) {
    if (this->specifier_.derivOrder == 0) {
      result_[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, center}] =
          Eigen::MatrixXd::Zero(this->dim1_ * this->dim1_, this->dim2_ * this->dim2_);
      resultPtr_.push_back(result_[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, center}].data());
    }
    else {
      for (auto const& derivKey :
           {Utils::Integrals::DerivKey::x, Utils::Integrals::DerivKey::y, Utils::Integrals::DerivKey::z}) {
        result_[{Utils::Integrals::Component::none, derivKey, center}] =
            Eigen::MatrixXd::Zero(this->dim1_ * this->dim1_, this->dim2_ * this->dim2_);
        resultPtr_.push_back(result_[{Utils::Integrals::Component::none, derivKey, center}].data());
      }
    }
  }

  if (!specifier.totalMass.has_value()) {
    throw std::runtime_error("No total Mass given in integral specifier.");
  }
  totalMass_ = specifier.totalMass.get();

  // Evaluate derivatives of overlap integrals:
  Utils::Integrals::IntegralSpecifier momentumSpecifier;
  momentumSpecifier.op = Utils::Integrals::Operator::Overlap;
  momentumSpecifier.derivOrder = 1;

  auto oneBodyInts = OneBodyIntegrals(scineBasis1, scineBasis1, momentumSpecifier);
  oneBodyInts.compute();
  vectorMomentumIntegrals_.emplace_back(oneBodyInts.getResult());
  if (scineBasis1 != scineBasis2) {
    auto oneBodyInts2 = OneBodyIntegrals(scineBasis2, scineBasis2, momentumSpecifier);
    oneBodyInts2.compute();
    vectorMomentumIntegrals_.emplace_back(oneBodyInts2.getResult());
  }
}

template class COMSaverDigester<IntegralSymmetry::twofold>;
template class COMSaverDigester<IntegralSymmetry::fourfold>;

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine
