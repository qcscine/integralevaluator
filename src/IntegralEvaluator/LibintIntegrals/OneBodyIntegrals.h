/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_ONEBODYINTEGRALS_H
#define INTEGRALEVALUATOR_ONEBODYINTEGRALS_H

/* internal */
#include <LibintIntegrals/LibintIntegrals.h>
/* external */
#include <Utils/DataStructures/IntegralSpecifier.h>
#include <boost/functional/hash.hpp>

namespace Scine {
namespace Integrals {

/**
 * @class OneBodyIntegrals
 * This class is called by the `LibintInterface`, and it computes the one-body integrals.
 */
class OneBodyIntegrals {
 private:
  IntegralEvaluatorMap result_;
  std::vector<Utils::Integrals::Component> relevantComponents_;
  std::vector<Utils::Integrals::DerivKey> relevantDerivKeys_;
  std::pair<size_t, size_t> dimension_;
  std::size_t numberOfOperators_;
  std::size_t numberOfResults_;
  std::size_t numberOfGeometricalDerivatives_;
  std::size_t numberOfCenters_;
  const Utils::Integrals::BasisSet& basis1_;
  const Utils::Integrals::BasisSet& basis2_;
  const Utils::Integrals::IntegralSpecifier& specifier_;

 public:
  /**
   * @brief Constructor for a one-body integral based on infos in `specifier`.
   * @param basis1
   * @param basis2
   * @param specifier
   */
  OneBodyIntegrals(const Utils::Integrals::BasisSet& basis1, const Utils::Integrals::BasisSet& basis2,
                   const Utils::Integrals::IntegralSpecifier& specifier);

  /**
   * @brief Do the actual computation
   */
  auto compute() -> void;

  /**
   * @brief After `compute` has been called, the result can be retrieved with this method.
   * @return
   */
  auto getResult() -> IntegralEvaluatorMap;

 private:
  auto _integral() -> void;

  auto _constructResultMap() -> void;
};

} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_ONEBODYINTEGRALS_H
