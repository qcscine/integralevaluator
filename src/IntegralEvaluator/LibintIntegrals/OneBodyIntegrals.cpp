/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <LibintIntegrals/BasisSetHandler.h>
#include <LibintIntegrals/Libint.h>
#include <LibintIntegrals/OneBodyIntegrals.h>

using namespace Scine;
using namespace Integrals;

OneBodyIntegrals::OneBodyIntegrals(const Utils::Integrals::BasisSet& basis1, const Utils::Integrals::BasisSet& basis2,
                                   const Utils::Integrals::IntegralSpecifier& specifier)
  : basis1_(basis1), basis2_(basis2), specifier_(specifier) {
  dimension_.first = basis1.nbf();
  dimension_.second = basis2.nbf();

  _constructResultMap();
}

auto OneBodyIntegrals::compute() -> void {
  _integral();
}

auto OneBodyIntegrals::_constructResultMap() -> void {
  libint2::Operator op = BasisSetHandler::scineToLibint(specifier_.op);

  switch (op) {
    case libint2::Operator::overlap:
      relevantComponents_ = {Utils::Integrals::Component::none};
      break;
    case libint2::Operator::kinetic:
      relevantComponents_ = {Utils::Integrals::Component::none};
      break;
    case libint2::Operator::nuclear:
      relevantComponents_ = {Utils::Integrals::Component::none};
      break;
    case libint2::Operator::emultipole1:
      relevantComponents_ = {Utils::Integrals::Component::x, Utils::Integrals::Component::y, Utils::Integrals::Component::z};
      break;
    default:
      throw std::runtime_error("Operator not available in the one-body integral routine!");
  }

  numberOfOperators_ = relevantComponents_.size();

  if (specifier_.derivOrder == 1) {
    relevantDerivKeys_ = {Utils::Integrals::DerivKey::x, Utils::Integrals::DerivKey::y, Utils::Integrals::DerivKey::z};
  }
  else {
    relevantDerivKeys_ = {
        Utils::Integrals::DerivKey::value,
    };
  }

  // The first center is the bra, the second center is the ket.
  // For the PointCharges operator, the remaining centers are the centers of the atoms (shells in the general case).
  if (specifier_.derivOrder > 0 && specifier_.op == Utils::Integrals::Operator::PointCharges) {
    assert(specifier_.atoms.has_value());
    numberOfCenters_ = specifier_.atoms.get().size() + 2;
    numberOfGeometricalDerivatives_ =
        libint2::num_geometrical_derivatives(specifier_.atoms.get().size(), specifier_.derivOrder);
  }
  else {
    numberOfCenters_ = (specifier_.derivOrder > 0) ? 2 : 1;
  }
  numberOfResults_ = numberOfOperators_ * numberOfCenters_ * relevantDerivKeys_.size();

  result_.reserve(numberOfResults_);

  for (auto const& component : relevantComponents_) {
    for (std::size_t center = 0; center < numberOfCenters_; ++center) {
      for (auto const& derivKey : relevantDerivKeys_) {
        result_[{component, derivKey, center}] = Eigen::MatrixXd::Zero(dimension_.first, dimension_.second);
      }
    }
  }
}

auto OneBodyIntegrals::_integral() -> void {
  double scaling = 1;

  libint2::Operator op = BasisSetHandler::scineToLibint(specifier_.op);

  auto s_engine = Libint::getEngine(basis1_, basis2_, op, specifier_.derivOrder);

  if (op == libint2::Operator::nuclear) {
    scaling = specifier_.typeVector[0].charge;
    // libint assumes negatively charged electrons with charge=(-1) and positive point charges.
    // therefore, we account for nuclear-nuclear repulsion with a ``-1``.
    // In the case of negatively charged particles, we have to take the absolute since the ``-`` is already implied.
    if (scaling > 0)
      scaling *= (-1);
    else
      scaling = std::abs(scaling);
    if (!specifier_.atoms.has_value())
      throw std::runtime_error("No atoms given in integral specifier.");
    auto libintAtoms = BasisSetHandler::scineToLibint(specifier_.atoms.get());
    auto pointCharges = libint2::make_point_charges(libintAtoms);
    s_engine.set_params(pointCharges);
  }
  else if (op == libint2::Operator::kinetic) {
    scaling = 1. / specifier_.typeVector[0].mass;
    if (specifier_.op == Utils::Integrals::Operator::KineticCOM) {
      if (!specifier_.totalMass.has_value())
        throw std::runtime_error("No total Mass given in integral specifier.");
      scaling -= 1 / specifier_.totalMass.get();
    }
  }
  else if (op == libint2::Operator::emultipole1) {
    scaling = specifier_.typeVector[0].charge;
    // libint assumes negatively charged electrons with charge=(-1) and positive point charges.
    // therefore, we account for nuclear-nuclear repulsion with a ``-1``.
    // In the case of negatively charged particles, we have to take the absolute since the ``-`` is already implied.
    if (scaling > 0)
      scaling *= (-1);
    else
      scaling = std::abs(scaling);
    if (!specifier_.multipoleOrigin.has_value())
      throw std::runtime_error("No multipole origin given in integral specifier.");
    s_engine.set_params(std::array<double, 3>{specifier_.multipoleOrigin.get()[0], specifier_.multipoleOrigin.get()[1],
                                              specifier_.multipoleOrigin.get()[2]});
  }

  const auto& buf_vec = s_engine.results(); // will point to computed shell sets --> const auto& is very important

  auto shell2bf1 = basis1_.shell2bf();
  auto shell2bf2 = basis2_.shell2bf();

  // This is the most delicate part: retrieving the correct indices.
  // iter1 and iter2 store the index of shell1 and shell2
  for (size_t s1 = 0; s1 < basis1_.size(); ++s1) {
    auto shell1 = BasisSetHandler::scineToLibint(basis1_[s1]);
    for (size_t s2 = 0; s2 < basis2_.size(); ++s2) {
      auto shell2 = BasisSetHandler::scineToLibint(basis2_[s2]);

      s_engine.compute(shell1, shell2);

      // Number of functions in shell 1 and 2. This depends on the angular momentum of the shell.
      auto n1 = shell1.size();
      auto bf1 = shell2bf1.at(s1);
      auto n2 = shell2.size();
      auto bf2 = shell2bf2.at(s2);

      for (auto const& component : relevantComponents_) {
        for (std::size_t center = 0; center < numberOfCenters_; ++center) {
          for (auto const& derivKey : relevantDerivKeys_) {
            auto index = static_cast<int>(component) * relevantDerivKeys_.size() * numberOfCenters_ +
                         center * relevantDerivKeys_.size() + static_cast<int>(derivKey);
            const auto* ints_shellset = buf_vec[index];
            // nullptr returned if the entire shell-set was screened out
            if (ints_shellset != nullptr) {
              Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>> tmp(buf_vec[index], n1, n2);
              result_[{component, derivKey, center}].block(bf1, bf2, n1, n2) = tmp * scaling;
              // for (std::size_t f1 = 0; f1 < n1; ++f1) {
              //  for (std::size_t f2 = 0; f2 < n2; ++f2) {
              //    // Eigen::Map<const Eigen::MatrixXd> buf_mat(ints_shellset, n1, n2);
              //    result_[{component, derivKey, center}](bf1 + f1, bf2 + f2) = ints_shellset[f1 * n2 + f2] * scaling;
              //  }
              //}
            }
          }
        }
      }
    } // s2
  }   // s1
}

auto OneBodyIntegrals::getResult() -> IntegralEvaluatorMap {
  return std::move(result_);
}
