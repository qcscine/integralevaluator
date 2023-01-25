/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal includes */
#include <LibintIntegrals/BasisSetHandler.h>
#include <LibintIntegrals/IntegralEvaluatorSettings.h>
#include <LibintIntegrals/Libint.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <LibintIntegrals/OneBodyIntegrals.h>
#include <LibintIntegrals/TwoBodyIntegrals/COMSaverDigester.h>
#include <LibintIntegrals/TwoBodyIntegrals/CauchySchwarzDensityPrescreener.h>
#include <LibintIntegrals/TwoBodyIntegrals/Evaluator.h>
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeDigester.h>
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombDigester.h>
#include <LibintIntegrals/TwoBodyIntegrals/SaverDigester.h>
/* External includes */
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Settings.h>

using namespace Scine;
using namespace Integrals;

LibintIntegrals::LibintIntegrals() {
  _settings = std::make_unique<IntegralEvaluatorSettings>();
}

LibintIntegrals::~LibintIntegrals() {
}

auto LibintIntegrals::initializeBasisSet(const std::string& name, const Utils::AtomCollection& atoms,
                                         bool generateShellPairs) -> Utils::Integrals::BasisSet {
  std::vector<libint2::Atom> libint_atoms = BasisSetHandler::scineToLibint(atoms);

  auto libint_basis = libint2::BasisSet(name, libint_atoms, true);

  libint_basis.set_pure(_settings->getBool("use_pure_spherical"));

  Utils::Integrals::BasisSet scineBasis =
      BasisSetHandler::libintToScine(libint_basis, atoms, _settings->getBool("use_pure_spherical"));

  if (generateShellPairs) {
    BasisSetHandler::generateShellPairs(scineBasis);
  }

  return scineBasis;
}

auto LibintIntegrals::initializeBasisSet(const std::unordered_map<Utils::ElementType, std::string>& names,
                                         const Utils::AtomCollection& atoms, bool generateShellPairs)
    -> Utils::Integrals::BasisSet {
  auto scineBasis = Utils::Integrals::BasisSet();

  for (auto const& elem_name_pair : names) {
    Utils::AtomCollection tmp_atom_coll;
    bool found = false;
    for (auto const& atom : atoms) {
      if (elem_name_pair.first == atom.getElementType()) {
        tmp_atom_coll.push_back(atom);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error("Element not contained in molecular structure");
    }
    std::vector<libint2::Atom> libint_atoms = BasisSetHandler::scineToLibint(tmp_atom_coll);
    auto libint_basis = libint2::BasisSet(elem_name_pair.second, libint_atoms, true);
    libint_basis.set_pure(_settings->getBool("use_pure_spherical"));
    scineBasis.append(BasisSetHandler::libintToScine(libint_basis, tmp_atom_coll, _settings->getBool("use_pure_spherical")));
  }

  if (generateShellPairs) {
    BasisSetHandler::generateShellPairs(scineBasis);
  }
  return scineBasis;
}

auto LibintIntegrals::evaluateOneBody(const Utils::Integrals::IntegralSpecifier& specifier,
                                      const Utils::Integrals::BasisSet& basis1, const Utils::Integrals::BasisSet& basis2)
    -> IntegralEvaluatorMap {
  auto oneBodyInts = OneBodyIntegrals(basis1, basis2, specifier);
  oneBodyInts.compute();
  return oneBodyInts.getResult();
}

auto LibintIntegrals::evaluateTwoBody(const Utils::Integrals::IntegralSpecifier& specifier,
                                      const Utils::Integrals::BasisSet& basis1, const Utils::Integrals::BasisSet& basis2)
    -> IntegralEvaluatorMap {
  if (!basis1.areShellPairsEvaluated()) {
    throw std::runtime_error("Evaluate shell pairs before performing the two-body integral evaluation!");
  }
  if (!basis2.areShellPairsEvaluated()) {
    throw std::runtime_error("Evaluate shell pairs before performing the two-body integral evaluation!");
  }

  auto prescreener = TwoBody::VoidPrescreener();

  IntegralEvaluatorMap result;

  if (specifier.op == Utils::Integrals::Operator::Coulomb) {
    if (basis1 == basis2) {
      auto saver = TwoBody::SaverDigester<TwoBody::IntegralSymmetry::eightfold>(basis1, basis2, specifier);
      auto eval = TwoBody::Evaluator<TwoBody::SaverDigester<TwoBody::IntegralSymmetry::eightfold>>(
          basis1, basis2, specifier, std::move(saver), std::move(prescreener));
      eval.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
      result = eval.getResult();
    }
    else {
      auto saver = TwoBody::SaverDigester<TwoBody::IntegralSymmetry::fourfold>(basis1, basis2, specifier);
      auto eval = TwoBody::Evaluator<TwoBody::SaverDigester<TwoBody::IntegralSymmetry::fourfold>>(
          basis1, basis2, specifier, std::move(saver), std::move(prescreener));
      eval.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
      result = eval.getResult();
    }
  }
  else if (specifier.op == Utils::Integrals::Operator::CoulombCOM) {
    if (basis1 == basis2) {
      auto saver = TwoBody::COMSaverDigester<TwoBody::IntegralSymmetry::fourfold>(basis1, basis2, specifier);
      auto eval = TwoBody::Evaluator<TwoBody::COMSaverDigester<TwoBody::IntegralSymmetry::fourfold>>(
          basis1, basis2, specifier, std::move(saver), std::move(prescreener));
      eval.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
      result = eval.getResult();
    }
    else {
      auto saver = TwoBody::COMSaverDigester<TwoBody::IntegralSymmetry::twofold>(basis1, basis2, specifier);
      auto eval = TwoBody::Evaluator<TwoBody::COMSaverDigester<TwoBody::IntegralSymmetry::twofold>>(
          basis1, basis2, specifier, std::move(saver), std::move(prescreener));
      eval.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
      result = eval.getResult();
    }
  }
  else {
    throw std::runtime_error("Invalid two-body operator");
  }

  return result;
}

std::string LibintIntegrals::name() {
  return std::string(model);
}

Utils::Settings& LibintIntegrals::settings() {
  return *_settings;
}

const Utils::Settings& LibintIntegrals::settings() const {
  return *_settings;
}

auto LibintIntegrals::evaluate(const Utils::Integrals::IntegralSpecifier& specifier, const Utils::Integrals::BasisSet& basis1,
                               const Utils::Integrals::BasisSet& basis2) -> IntegralEvaluatorMap {
  if (specifier.op == Utils::Integrals::Operator::Kinetic || specifier.op == Utils::Integrals::Operator::KineticCOM ||
      specifier.op == Utils::Integrals::Operator::PointCharges || specifier.op == Utils::Integrals::Operator::Overlap ||
      specifier.op == Utils::Integrals::Operator::Dipole) {
    return evaluateOneBody(specifier, basis1, basis2);
  }
  else if (specifier.op == Utils::Integrals::Operator::Coulomb || specifier.op == Utils::Integrals::Operator::CoulombCOM) {
    return evaluateTwoBody(specifier, basis1, basis2);
  }
  else
    throw std::runtime_error("Operator not implemented!");
}
auto LibintIntegrals::generateShellPairs(Utils::Integrals::BasisSet& basis, bool performOverlapPrescreening,
                                         double threshold, bool calculateCauchySchwarzFactor) -> void {
  if (basis.areShellPairsEvaluated()) {
    std::cout << "Shell pairs already evaluated. Nothing was done" << std::endl;
    return;
  }
  BasisSetHandler::generateShellPairs(basis, performOverlapPrescreening, threshold, calculateCauchySchwarzFactor);
}

auto LibintIntegrals::evaluateTwoBodyDirectBo(const Utils::Integrals::IntegralSpecifier& specifier,
                                              const Utils::Integrals::BasisSet& basis1, const Utils::Integrals::BasisSet& basis2,
                                              const Utils::DensityMatrix& dm1, double prescreeningThreshold)
    -> std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix> {
  if (basis1 == basis2) {
    auto prescreener = Integrals::TwoBody::CauchySchwarzDensityPrescreener(basis1, dm1, prescreeningThreshold);

    auto saver = Integrals::TwoBody::CoulombExchangeDigester(basis1, basis2, specifier, dm1);
    auto evaluator =
        Integrals::TwoBody::Evaluator<Integrals::TwoBody::CoulombExchangeDigester, Integrals::TwoBody::CauchySchwarzDensityPrescreener>(
            basis1, basis2, specifier, std::move(saver), std::move(prescreener));

    evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();

    Utils::SpinAdaptedMatrix J = evaluator.getResult().first;
    Utils::SpinAdaptedMatrix K = evaluator.getResult().second;

    return {J, K};
  }
  else {
    auto prescreener = Integrals::TwoBody::VoidPrescreener();

    auto saver = Integrals::TwoBody::CoulombExchangeDigester(basis1, basis2, specifier, dm1);
    auto evaluator =
        Integrals::TwoBody::Evaluator<Integrals::TwoBody::CoulombExchangeDigester, Integrals::TwoBody::VoidPrescreener>(
            basis1, basis2, specifier, std::move(saver), std::move(prescreener));

    evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();

    Utils::SpinAdaptedMatrix J = evaluator.getResult().first;
    Utils::SpinAdaptedMatrix K = evaluator.getResult().second;

    return {J, K};
  }
}

auto LibintIntegrals::evaluateTwoBodyDirectPreBo(const Utils::Integrals::IntegralSpecifier& specifier,
                                                 const Utils::Integrals::BasisSet& basis1,
                                                 const Utils::Integrals::BasisSet& basis2,
                                                 const Utils::DensityMatrix& dm1, const Utils::DensityMatrix& dm2)
    -> std::pair<Eigen::MatrixXd, Eigen::MatrixXd> {
  if (specifier.typeVector.at(0).symbol == specifier.typeVector.at(1).symbol) {
    throw std::runtime_error("BO routine was called with different particle types.");
  }

  auto prescreener = Integrals::TwoBody::VoidPrescreener();

  auto saver = Integrals::TwoBody::TwoTypeCoulombDigester(basis1, basis2, specifier, dm1, dm2);
  auto evaluator = Integrals::TwoBody::Evaluator<Integrals::TwoBody::TwoTypeCoulombDigester>(
      basis1, basis2, specifier, std::move(saver), std::move(prescreener));
  evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();

  return {evaluator.getResult().first, evaluator.getResult().second};
}
