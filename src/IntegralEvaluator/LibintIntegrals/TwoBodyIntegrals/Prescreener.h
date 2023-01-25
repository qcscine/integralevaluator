/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INTEGRALEVALUATOR_PRESCREENER_H
#define INTEGRALEVALUATOR_PRESCREENER_H

namespace Scine {

namespace Utils {
namespace Integrals {
class BasisSet;
}
} // namespace Utils

namespace Utils {
class DensityMatrix;
} // namespace Utils

namespace Integrals {
namespace TwoBody {
/**
 * @class TwoBodiesIntegralsPrescreener @file TwoBodiesIntegralsPrescreener
 * @brief CRTP functor representing a static interface for a prescreener.
 * @tparam CrtpClass The CRTP parameter.
 */
template<typename CrtpClass>
class TwoBodiesIntegralsPrescreener {
 public:
  explicit TwoBodiesIntegralsPrescreener(const Utils::Integrals::BasisSet& scineBasis, double prescreenThreshold = 1e-12);

  /**
   * @brief Function to be called in order to use a prescreener.
   * @returns A bool flagging the validity of the shell set.
   */
  bool operator()(int shell1, int shell2, int shell3, int shell4, double shellSetCriterion) const {
    return derived().isSignificantImpl(shell1, shell2, shell3, shell4, shellSetCriterion);
  }

 protected:
  const Utils::Integrals::BasisSet& scineBasis_;
  const double prescreeningThreshold_;

 private:
  CrtpClass& derived() {
    return static_cast<CrtpClass&>(*this);
  }
  const CrtpClass& derived() const {
    return static_cast<const CrtpClass&>(*this);
  }
};

template<typename CrtpClass>
TwoBodiesIntegralsPrescreener<CrtpClass>::TwoBodiesIntegralsPrescreener(const Utils::Integrals::BasisSet& scineBasis,
                                                                        double prescreenThreshold)
  : scineBasis_(scineBasis), prescreeningThreshold_(prescreenThreshold) {
}

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_PRESCREENER_H
