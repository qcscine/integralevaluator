/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_INTEGRALEVALUATORSETTINGS_H
#define INTEGRALEVALUATOR_INTEGRALEVALUATORSETTINGS_H

#include <Utils/Settings.h>

namespace Scine {
namespace Integrals {

struct IntegralEvaluatorSettings : public Utils::Settings {
  IntegralEvaluatorSettings() : Utils::Settings("Settings for the Integral Evaluator.") {
    Utils::UniversalSettings::BoolDescriptor pure_spherical(
        "If true, use only solid harmonics (i.e., spherical harmonics), else, use only cartesians. Mixing is not "
        "possible.");
    pure_spherical.setDefaultValue(true);
    _fields.push_back("use_pure_spherical", std::move(pure_spherical));

    resetToDefaults();
  }
};

} // namespace Integrals
} // namespace Scine
#endif // INTEGRALEVALUATOR_INTEGRALEVALUATORSETTINGS_H
