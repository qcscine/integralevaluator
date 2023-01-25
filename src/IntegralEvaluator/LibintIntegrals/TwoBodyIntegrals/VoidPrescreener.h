/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_VOIDPRESCREENER_H
#define INTEGRALEVALUATOR_VOIDPRESCREENER_H

#include <LibintIntegrals/TwoBodyIntegrals/Prescreener.h>
#include <Utils/DataStructures/BasisSet.h>

namespace Scine {
namespace Integrals {
namespace TwoBody {

class VoidPrescreener : public TwoBodiesIntegralsPrescreener<VoidPrescreener> {
 public:
  VoidPrescreener() : TwoBodiesIntegralsPrescreener<VoidPrescreener>(voidScineBasis_){};
  VoidPrescreener(VoidPrescreener&&) = default;
  VoidPrescreener& operator=(VoidPrescreener&&) = delete;
  VoidPrescreener(const VoidPrescreener&) = delete;
  VoidPrescreener& operator=(const VoidPrescreener&) = delete;

  bool isSignificantImpl(int, int, int, int, double) const {
    return true;
  }

 private:
  Utils::Integrals::BasisSet voidScineBasis_;
};

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_VOIDPRESCREENER_H
