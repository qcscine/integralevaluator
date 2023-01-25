/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <LibintIntegrals/TwoBodyIntegrals/SaverDigester.h>
#include <LibintIntegrals/TwoBodyIntegrals/SymmetryHelper.h>

#define UNUSED(expr) \
  do {               \
    (void)(expr);    \
  } while (0)

namespace Scine {
namespace Integrals {
namespace TwoBody {

template<>
void SaverDigester<IntegralSymmetry::eightfold>::digestImpl(double integralValue, int i, int j, int k, int l, int index,
                                                            double degeneracy) {
  UNUSED(degeneracy);

  auto& ptr_data = resultPtr_[index];

  integralValue *= this->scaling_;

  auto ij = i * dim1_ + j;
  auto ji = j * dim1_ + i;
  auto kl = k * dim2_ + l;
  auto lk = l * dim2_ + k;

  // 4-fold symmetry
  // ij,kl = ij,lk = ji,kl = ji,lk
  ptr_data[ij + dim1sq_ * (kl)] = integralValue;
  ptr_data[ij + dim1sq_ * (lk)] = integralValue;
  ptr_data[ji + dim1sq_ * (kl)] = integralValue;
  ptr_data[ji + dim1sq_ * (lk)] = integralValue;
  // 8-fold symmetry:
  // kl,ij  = lk,ji = lk,ij  = kl,ji
  ptr_data[kl + dim1sq_ * (ij)] = integralValue;
  ptr_data[kl + dim1sq_ * (ji)] = integralValue;
  ptr_data[lk + dim1sq_ * (ij)] = integralValue;
  ptr_data[lk + dim1sq_ * (ji)] = integralValue;
}

template<>
void SaverDigester<IntegralSymmetry::fourfold>::digestImpl(double integralValue, int i, int j, int k, int l, int index,
                                                           double degeneracy) {
  UNUSED(degeneracy);

  auto& ptr_data = resultPtr_[index];

  integralValue *= this->scaling_;

  // 4-fold symmetry
  // ij,kl = ij,lk = ji,kl = ji,lk
  ptr_data[i * dim1_ + j + dim1sq_ * (k * dim2_ + l)] = integralValue;
  ptr_data[i * dim1_ + j + dim1sq_ * (l * dim2_ + k)] = integralValue;
  ptr_data[j * dim1_ + i + dim1sq_ * (k * dim2_ + l)] = integralValue;
  ptr_data[j * dim1_ + i + dim1sq_ * (l * dim2_ + k)] = integralValue;
}

template<IntegralSymmetry symmetry>
const IntegralEvaluatorMap& SaverDigester<symmetry>::getResultImpl() const {
  return result_;
}

template<IntegralSymmetry symmetry>
double SaverDigester<symmetry>::computeDegeneracyImpl(int shell1, int shell2, int shell3, int shell4) {
  UNUSED(shell1);
  UNUSED(shell2);
  UNUSED(shell3);
  UNUSED(shell4);
  return 1;
}

template<IntegralSymmetry symmetry>
void SaverDigester<symmetry>::initializeImpl(int numberThreads) {
  // Parallelizing the SaverDigester will not be beneficial.
  UNUSED(numberThreads);
}

template<IntegralSymmetry symmetry>
void SaverDigester<symmetry>::finalizeImpl() {
}

template<IntegralSymmetry symmetry>
SaverDigester<symmetry>::SaverDigester(const Utils::Integrals::BasisSet& scineBasis1,
                                       const Utils::Integrals::BasisSet& scineBasis2,
                                       const Utils::Integrals::IntegralSpecifier& specifier)
  : Digester<SaverDigester<symmetry>>(scineBasis1, scineBasis2, specifier) {
  static_assert(symmetry == IntegralSymmetry::eightfold || symmetry == IntegralSymmetry::fourfold,
                "Symmetry must be either 4- or 8-fold.");

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
}

template class SaverDigester<IntegralSymmetry::fourfold>;
template class SaverDigester<IntegralSymmetry::eightfold>;

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine
