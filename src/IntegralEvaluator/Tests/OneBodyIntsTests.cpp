/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <LibintIntegrals/BasisSetHandler.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <Utils/Constants.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>
#include <Eigen/Eigenvalues>

// Only for reference test
#include <Utils/Constants.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/DensityMatrixBuilder.h>

using namespace Scine;
using namespace Integrals;

using namespace testing;

class OneBodyIntsTest : public Test {};

TEST_F(OneBodyIntsTest, TestOverlap) {
  // Reference
  Eigen::MatrixXd pyScfOverlapH2Def2SVP;
  pyScfOverlapH2Def2SVP.resize(10, 10);
  pyScfOverlapH2Def2SVP << 1.00000000e+00, 6.84799825e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      2.53303682e-01, 4.14114249e-01, -2.86782459e-01, 0.00000000e+00, 0.00000000e+00, 6.84799825e-01, 1.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 4.14114249e-01, 7.30845790e-01, -1.73676647e-01, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.86782459e-01,
      1.73676647e-01, -3.98093618e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.27845428e-01, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 1.27845428e-01, 2.53303682e-01, 4.14114249e-01, 2.86782459e-01, 0.00000000e+00,
      0.00000000e+00, 1.00000000e+00, 6.84799825e-01, -8.49505555e-18, 0.00000000e+00, 0.00000000e+00, 4.14114249e-01,
      7.30845790e-01, 1.73676647e-01, 0.00000000e+00, 0.00000000e+00, 6.84799825e-01, 1.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, -2.86782459e-01, -1.73676647e-01, -3.98093618e-01, 0.00000000e+00, 0.00000000e+00,
      -8.49505555e-18, 0.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 1.27845428e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.27845428e-01, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00;

  std::stringstream h2("2\n\n"
                       "H 0 0 0\n"
                       "H 1.2 0 0");
  auto scineAtoms = Utils::XyzStreamHandler::read(h2);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis = eval.initializeBasisSet("def2-svp", scineAtoms);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Overlap;

  auto result_map = LibintIntegrals::evaluate(specifier, basis, basis);

  const auto& result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  for (int row = 0; row < 10; ++row) {
    for (int col = 0; col < 10; ++col) {
      EXPECT_THAT(result(row, col), DoubleNear(pyScfOverlapH2Def2SVP(row, col), 1e-8));
    }
  }
}

TEST_F(OneBodyIntsTest, TestOverlap2) {
  std::stringstream h2os("4\n\n"
                         "H 0 0 0\n"
                         "H 1.2 0 0\n"
                         "O -1 0 0\n"
                         "S 2 0 0");

  auto scineAtoms = Utils::XyzStreamHandler::read(h2os);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis = eval.initializeBasisSet(
      {{Utils::ElementType::H, "def2-svp"}, {Utils::ElementType::O, "def2-svp"}, {Utils::ElementType::S, "def2-svp"}},
      scineAtoms);

  auto basis_ref = eval.initializeBasisSet("def2-svp", scineAtoms);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Overlap;

  auto result_map = LibintIntegrals::evaluate(specifier, basis, basis);

  auto result_map_ref = LibintIntegrals::evaluate(specifier, basis_ref, basis_ref);

  const auto& result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  const auto& result_ref = result_map_ref[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  EXPECT_EQ(result.rows(), result_ref.rows());
  EXPECT_EQ(result.cols(), result_ref.cols());

  auto shell2bf = basis.shell2bf();
  auto shell2bf_ref = basis_ref.shell2bf();

  for (size_t s1 = 0; s1 < basis.size(); ++s1) {
    for (size_t s2 = 0; s2 < basis_ref.size(); ++s2) {
      if (basis[s1].getVecAlpha() == basis_ref[s2].getVecAlpha()) {
        auto n1 = basis[s1].size();
        auto bf1 = shell2bf.at(s1);
        auto bf2 = shell2bf_ref.at(s2);
        auto n2 = basis_ref[s2].size();
        assert(n1 == n2);
        for (std::size_t f1 = 0; f1 < n1; ++f1) {
          for (std::size_t f2 = 0; f2 < n2; ++f2) {
            EXPECT_THAT(result(bf1 + f1, bf1 + f2), DoubleNear(result_ref(bf2 + f1, bf2 + f2), 1e-8));
          }
        }
      }
    }
  }
}

TEST_F(OneBodyIntsTest, TestOverlap3) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.0 0.0 0.0\n"
                             "H    0.740848 0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis = eval.initializeBasisSet("cc-pvtz", atoms.first);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Overlap;

  auto result_map = LibintIntegrals::evaluate(specifier, basis, basis);

  const auto& result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  for (int row = 0; row < result.rows(); ++row) {
    for (int col = 0; col < row; ++col) {
      EXPECT_THAT(result(row, col), DoubleNear(result(col, row), 1e-12));
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(result.selfadjointView<Eigen::Lower>());

  for (int row = 0; row < es.eigenvalues().size(); ++row) {
    EXPECT_TRUE(es.eigenvalues()[row] > 0);
  }
}

TEST_F(OneBodyIntsTest, TestOverlapDerivative) {
  // Reference

  std::stringstream h2("2\n\n"
                       "H 0 0 0\n"
                       "H 1.2 0 0");
  auto atoms = Utils::XyzStreamHandler::read(h2);

  // Central difference:
  std::stringstream h2_h0("2\n\n"
                          "H 0 0 0\n"
                          "H 1.2001 0 0");
  auto atoms_h0 = Utils::XyzStreamHandler::read(h2_h0);

  std::stringstream h2_h1("2\n\n"
                          "H 0 0 0\n"
                          "H 1.1999 0 0");
  auto atoms_h1 = Utils::XyzStreamHandler::read(h2_h1);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  std::string name = "def2-svp";

  auto basis = eval.initializeBasisSet(name, atoms);
  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Overlap;
  specifier.derivOrder = 1;

  auto result_map = LibintIntegrals::evaluate(specifier, basis, basis);
  auto result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 0}];
  result += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 1}];

  Utils::Integrals::IntegralSpecifier testSpecifier;
  testSpecifier.op = Utils::Integrals::Operator::Overlap;
  std::vector<Eigen::MatrixXd> centralDifference;

  // Derivative in x-direction:
  auto tmpBasis = eval.initializeBasisSet(name, atoms_h0);
  auto tmp_map = LibintIntegrals::evaluate(testSpecifier, tmpBasis, tmpBasis);
  auto tmp_res = tmp_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  centralDifference.push_back(tmp_res);
  tmpBasis = eval.initializeBasisSet(name, atoms_h1);
  tmp_map = LibintIntegrals::evaluate(testSpecifier, tmpBasis, tmpBasis);
  tmp_res = tmp_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  centralDifference.push_back(tmp_res);

  Eigen::MatrixXd xDerivative =
      1 / (0.0002 * Utils::Constants::bohr_per_angstrom) * (centralDifference[0] - centralDifference[1]);

  auto shell2bf = basis.shell2bf();
  for (size_t s1 = 0; s1 < basis.size(); ++s1) {
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      auto n1 = basis[s1].size();
      auto bf1 = shell2bf.at(s1);
      auto bf2 = shell2bf.at(s2);
      auto n2 = basis[s2].size();

      for (std::size_t f1 = 0; f1 < n1; ++f1) {
        for (std::size_t f2 = 0; f2 < n2; ++f2) {
          if (basis[s1].getShift() == basis[s2].getShift()) {
            auto tmp = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 0}](bf1 + f1, bf2 + f2);
            tmp += result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 1}](bf1 + f1, bf2 + f2);
            EXPECT_THAT(xDerivative(bf1 + f1, bf2 + f2), DoubleNear(tmp, 1e-6));
          }
          else {
            auto tmp = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::x, 0}](bf1 + f1, bf2 + f2);
            EXPECT_THAT(xDerivative(bf1 + f1, bf2 + f2), DoubleNear(tmp, 1e-6));
          }
        }
      }
    } // s2
  }   // s1
}

TEST_F(OneBodyIntsTest, TestCore) {
  // Reference
  Eigen::MatrixXd pyScfHCoreH2Def2SVP;
  pyScfHCoreH2Def2SVP.resize(10, 10);
  pyScfHCoreH2Def2SVP << -0.84872143, -0.74618239, -0.11906097, 0., 0., -0.39465425, -0.51067364, 0.34020547, 0., 0.,
      -0.74618239, -0.76538332, -0.10274864, 0., 0., -0.51067364, -0.63111362, 0.21472884, 0., 0., -0.11906097,
      -0.10274864, 0.55434149, 0., 0., -0.34020547, -0.21472884, 0.07485182, 0., 0., 0., 0., 0., 0.63428204, 0., 0., 0.,
      0., -0.12945872, 0., 0., 0., 0., 0., 0.63428204, 0., 0., 0., 0., -0.12945872, -0.39465425, -0.51067364,
      -0.34020547, 0., 0., -0.84872143, -0.74618239, 0.11906097, 0., 0., -0.51067364, -0.63111362, -0.21472884, 0., 0.,
      -0.74618239, -0.76538332, 0.10274864, 0., 0., 0.34020547, 0.21472884, 0.07485182, 0., 0., 0.11906097, 0.10274864,
      0.55434149, 0., 0., 0., 0., 0., -0.12945872, 0., 0., 0., 0., 0.63428204, 0., 0., 0., 0., 0., -0.12945872, 0., 0.,
      0., 0., 0.63428204;

  std::stringstream h2("2\n\n"
                       "H  0.0 0.0 0.0\n"
                       "H  1.2 0.0 0.0");

  auto scineAtoms = Utils::XyzStreamHandler::read(h2);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis = eval.initializeBasisSet("def2-svp", scineAtoms);

  Utils::Integrals::IntegralSpecifier kinetic;
  kinetic.op = Utils::Integrals::Operator::Kinetic;

  auto kinetic_result_map = LibintIntegrals::evaluate(kinetic, basis, basis);

  Utils::Integrals::IntegralSpecifier point_charges;
  point_charges.op = Utils::Integrals::Operator::PointCharges;
  point_charges.atoms = scineAtoms;

  auto point_charges_result_map = LibintIntegrals::evaluate(point_charges, basis, basis);

  Eigen::MatrixXd result = kinetic_result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  result += point_charges_result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  for (int row = 0; row < 10; ++row) {
    for (int col = 0; col < 10; ++col) {
      EXPECT_THAT(result(row, col), DoubleNear(pyScfHCoreH2Def2SVP(row, col), 1e-8));
    }
  }

  // Utils::Integrals::IntegralSpecifier point_charges_deriv = point_charges;
  // point_charges_deriv.derivOrder = 1;

  // auto point_charges_deriv_result_map = LibintIntegrals::evaluate(point_charges_deriv, basis, basis);
}
