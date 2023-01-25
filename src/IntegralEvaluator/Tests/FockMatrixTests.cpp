/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <LibintIntegrals/Libint.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <LibintIntegrals/TwoBodyIntegrals/CauchySchwarzDensityPrescreener.h>
#include <LibintIntegrals/TwoBodyIntegrals/Evaluator.h>
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeDigester.h>
#include <LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombDigester.h>
#include <Utils/Constants.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Scf/LcaoUtils/DensityMatrixBuilder.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>
#include <Eigen/Eigenvalues>

using namespace Scine;
using namespace Integrals;

using namespace testing;

class FockMatrixTest : public Test {};

TEST_F(FockMatrixTest, TestRestrictedJK) {
  omp_set_num_threads(1);

  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0\n"
                             "H    0.0 0.0 0.0");

  auto scineAtoms = Utils::XyzStreamHandler::read(xyzInput);

  // Reference
  Eigen::MatrixXd D;
  D.resize(10, 10);
  Eigen::MatrixXd J;
  J.resize(10, 10);
  Eigen::MatrixXd K;
  K.resize(10, 10);

  D << 0.5104938957, 0.0465016349, 0.0000000000, -0.0000000000, -0.0473914727, 0.5104938957, 0.0465016349, 0.0000000000,
      0.0000000000, 0.0473914727, 0.0465016349, 0.0042359019, 0.0000000000, -0.0000000000, -0.0043169585, 0.0465016349,
      0.0042359019, 0.0000000000, 0.0000000000, 0.0043169585, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000, -0.0000000000, -0.0000000000,
      -0.0000000000, -0.0473914727, -0.0043169585, -0.0000000000, 0.0000000000, 0.0043995662, -0.0473914727,
      -0.0043169585, -0.0000000000, -0.0000000000, -0.0043995662, 0.5104938957, 0.0465016349, 0.0000000000,
      -0.0000000000, -0.0473914727, 0.5104938957, 0.0465016349, 0.0000000000, 0.0000000000, 0.0473914727, 0.0465016349,
      0.0042359019, 0.0000000000, -0.0000000000, -0.0043169585, 0.0465016349, 0.0042359019, 0.0000000000, 0.0000000000,
      0.0043169585, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0473914727, 0.0043169585, 0.0000000000,
      -0.0000000000, -0.0043995662, 0.0473914727, 0.0043169585, 0.0000000000, 0.0000000000, 0.0043995662;

  J << 1.5492381470, 0.9095145472, 0.0000000000, -0.0000000000, -0.2116646284, 0.9723322859, 0.7892327421, 0.0000000000,
      0.0000000000, 0.9294739085, 0.9095145472, 0.9647712385, 0.0000000000, -0.0000000000, -0.1463753952, 0.7892327421,
      0.8910249847, 0.0000000000, 0.0000000000, 0.3109522093, 0.0000000000, 0.0000000000, 1.4101623569, -0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, 0.7496185965, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000,
      -0.0000000000, 1.4101623569, -0.0000000000, -0.0000000000, -0.0000000000, 0.0000000000, 0.7496185965,
      -0.0000000000, -0.2116646284, -0.1463753952, -0.0000000000, -0.0000000000, 1.4917172321, -0.9294739085,
      -0.3109522093, -0.0000000000, -0.0000000000, -0.4490206507, 0.9723322859, 0.7892327421, 0.0000000000,
      -0.0000000000, -0.9294739085, 1.5492381470, 0.9095145472, 0.0000000000, 0.0000000000, 0.2116646284, 0.7892327421,
      0.8910249847, 0.0000000000, -0.0000000000, -0.3109522093, 0.9095145472, 0.9647712385, 0.0000000000, 0.0000000000,
      0.1463753952, 0.0000000000, 0.0000000000, 0.7496185965, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      1.4101623569, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.7496185965, -0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 1.4101623569, 0.0000000000, 0.9294739085, 0.3109522093, 0.0000000000,
      -0.0000000000, -0.4490206507, 0.2116646284, 0.1463753952, 0.0000000000, 0.0000000000, 1.4917172321;

  K << 1.3285785199, 0.8753108146, 0.0000000000, -0.0000000000, -0.5225822918, 1.1833149770, 0.8494150509, 0.0000000000,
      0.0000000000, 0.6973039666, 0.8753108146, 0.6203480046, 0.0000000000, -0.0000000000, -0.3766166105, 0.8494150509,
      0.6144974143, 0.0000000000, 0.0000000000, 0.4101652644, 0.0000000000, 0.0000000000, 0.2352726804, -0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, 0.1764948067, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000,
      -0.0000000000, 0.2352726804, -0.0000000000, -0.0000000000, -0.0000000000, -0.0000000000, 0.1764948067,
      -0.0000000000, -0.5225822918, -0.3766166105, -0.0000000000, -0.0000000000, 0.4948970158, -0.6973039666,
      -0.4101652644, -0.0000000000, -0.0000000000, -0.2743137935, 1.1833149770, 0.8494150509, 0.0000000000,
      -0.0000000000, -0.6973039666, 1.3285785199, 0.8753108146, 0.0000000000, 0.0000000000, 0.5225822918, 0.8494150509,
      0.6144974143, 0.0000000000, -0.0000000000, -0.4101652644, 0.8753108146, 0.6203480046, 0.0000000000, 0.0000000000,
      0.3766166105, 0.0000000000, 0.0000000000, 0.1764948067, -0.0000000000, -0.0000000000, 0.0000000000, 0.0000000000,
      0.2352726804, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.1764948067, -0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 0.2352726804, 0.0000000000, 0.6973039666, 0.4101652644, 0.0000000000,
      -0.0000000000, -0.2743137935, 0.5225822918, 0.3766166105, 0.0000000000, 0.0000000000, 0.4948970158;

  auto coeffs = Eigen::MatrixXd::Random(10, 10);

  Utils::MolecularOrbitals mos = Utils::MolecularOrbitals::createFromRestrictedCoefficients(coeffs);

  Utils::LcaoUtils::DensityMatrixBuilder builder(mos);
  auto densityMatrix = builder.generateRestrictedForNumberElectrons(2);

  densityMatrix.setDensity(std::move(D), 2);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", true);

  auto basis = eval.initializeBasisSet("def2-svp", scineAtoms);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;

  auto prescreener = TwoBody::VoidPrescreener();

  std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix> result;

  auto saver = TwoBody::CoulombExchangeDigester(basis, basis, specifier, densityMatrix);
  auto evaluator = TwoBody::Evaluator<TwoBody::CoulombExchangeDigester>(basis, basis, specifier, std::move(saver),
                                                                        std::move(prescreener));
  evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
  result = evaluator.getResult();

  for (int row = 0; row < 10; ++row) {
    for (int col = 0; col < 10; ++col) {
      EXPECT_THAT(result.first.restrictedMatrix()(row, col), DoubleNear(J(row, col), 1e-8));
      EXPECT_THAT(result.second.restrictedMatrix()(row, col), DoubleNear(K(row, col), 1e-8));
    }
  }
}

TEST_F(FockMatrixTest, TestRestrictedJK_CauchySchwatz) {
  omp_set_num_threads(1);

  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0\n"
                             "H    0.0 0.0 0.0");

  auto scineAtoms = Utils::XyzStreamHandler::read(xyzInput);

  // Reference
  Eigen::MatrixXd D;
  D.resize(10, 10);
  Eigen::MatrixXd J;
  J.resize(10, 10);
  Eigen::MatrixXd K;
  K.resize(10, 10);

  D << 0.5104938957, 0.0465016349, 0.0000000000, -0.0000000000, -0.0473914727, 0.5104938957, 0.0465016349, 0.0000000000,
      0.0000000000, 0.0473914727, 0.0465016349, 0.0042359019, 0.0000000000, -0.0000000000, -0.0043169585, 0.0465016349,
      0.0042359019, 0.0000000000, 0.0000000000, 0.0043169585, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000, -0.0000000000, -0.0000000000,
      -0.0000000000, -0.0473914727, -0.0043169585, -0.0000000000, 0.0000000000, 0.0043995662, -0.0473914727,
      -0.0043169585, -0.0000000000, -0.0000000000, -0.0043995662, 0.5104938957, 0.0465016349, 0.0000000000,
      -0.0000000000, -0.0473914727, 0.5104938957, 0.0465016349, 0.0000000000, 0.0000000000, 0.0473914727, 0.0465016349,
      0.0042359019, 0.0000000000, -0.0000000000, -0.0043169585, 0.0465016349, 0.0042359019, 0.0000000000, 0.0000000000,
      0.0043169585, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0473914727, 0.0043169585, 0.0000000000,
      -0.0000000000, -0.0043995662, 0.0473914727, 0.0043169585, 0.0000000000, 0.0000000000, 0.0043995662;

  J << 1.5492381470, 0.9095145472, 0.0000000000, -0.0000000000, -0.2116646284, 0.9723322859, 0.7892327421, 0.0000000000,
      0.0000000000, 0.9294739085, 0.9095145472, 0.9647712385, 0.0000000000, -0.0000000000, -0.1463753952, 0.7892327421,
      0.8910249847, 0.0000000000, 0.0000000000, 0.3109522093, 0.0000000000, 0.0000000000, 1.4101623569, -0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, 0.7496185965, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000,
      -0.0000000000, 1.4101623569, -0.0000000000, -0.0000000000, -0.0000000000, 0.0000000000, 0.7496185965,
      -0.0000000000, -0.2116646284, -0.1463753952, -0.0000000000, -0.0000000000, 1.4917172321, -0.9294739085,
      -0.3109522093, -0.0000000000, -0.0000000000, -0.4490206507, 0.9723322859, 0.7892327421, 0.0000000000,
      -0.0000000000, -0.9294739085, 1.5492381470, 0.9095145472, 0.0000000000, 0.0000000000, 0.2116646284, 0.7892327421,
      0.8910249847, 0.0000000000, -0.0000000000, -0.3109522093, 0.9095145472, 0.9647712385, 0.0000000000, 0.0000000000,
      0.1463753952, 0.0000000000, 0.0000000000, 0.7496185965, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      1.4101623569, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.7496185965, -0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 1.4101623569, 0.0000000000, 0.9294739085, 0.3109522093, 0.0000000000,
      -0.0000000000, -0.4490206507, 0.2116646284, 0.1463753952, 0.0000000000, 0.0000000000, 1.4917172321;

  K << 1.3285785199, 0.8753108146, 0.0000000000, -0.0000000000, -0.5225822918, 1.1833149770, 0.8494150509, 0.0000000000,
      0.0000000000, 0.6973039666, 0.8753108146, 0.6203480046, 0.0000000000, -0.0000000000, -0.3766166105, 0.8494150509,
      0.6144974143, 0.0000000000, 0.0000000000, 0.4101652644, 0.0000000000, 0.0000000000, 0.2352726804, -0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, 0.1764948067, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000,
      -0.0000000000, 0.2352726804, -0.0000000000, -0.0000000000, -0.0000000000, -0.0000000000, 0.1764948067,
      -0.0000000000, -0.5225822918, -0.3766166105, -0.0000000000, -0.0000000000, 0.4948970158, -0.6973039666,
      -0.4101652644, -0.0000000000, -0.0000000000, -0.2743137935, 1.1833149770, 0.8494150509, 0.0000000000,
      -0.0000000000, -0.6973039666, 1.3285785199, 0.8753108146, 0.0000000000, 0.0000000000, 0.5225822918, 0.8494150509,
      0.6144974143, 0.0000000000, -0.0000000000, -0.4101652644, 0.8753108146, 0.6203480046, 0.0000000000, 0.0000000000,
      0.3766166105, 0.0000000000, 0.0000000000, 0.1764948067, -0.0000000000, -0.0000000000, 0.0000000000, 0.0000000000,
      0.2352726804, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.1764948067, -0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 0.2352726804, 0.0000000000, 0.6973039666, 0.4101652644, 0.0000000000,
      -0.0000000000, -0.2743137935, 0.5225822918, 0.3766166105, 0.0000000000, 0.0000000000, 0.4948970158;

  auto coeffs = Eigen::MatrixXd::Random(10, 10);
  Utils::MolecularOrbitals mos = Utils::MolecularOrbitals::createFromRestrictedCoefficients(coeffs);
  Utils::LcaoUtils::DensityMatrixBuilder builder(mos);
  auto densityMatrix = builder.generateRestrictedForNumberElectrons(2);
  densityMatrix.setDensity(std::move(D), densityMatrix.numberElectrons());

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", true);

  auto basis = eval.initializeBasisSet("def2-svp", scineAtoms);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;

  auto prescreener = TwoBody::CauchySchwarzDensityPrescreener(basis, densityMatrix);

  std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix> result;

  auto saver = TwoBody::CoulombExchangeDigester(basis, basis, specifier, densityMatrix);
  auto evaluator = TwoBody::Evaluator<TwoBody::CoulombExchangeDigester, TwoBody::CauchySchwarzDensityPrescreener>(
      basis, basis, specifier, std::move(saver), std::move(prescreener));
  evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
  result = evaluator.getResult();

  for (int row = 0; row < 10; ++row) {
    for (int col = 0; col < 10; ++col) {
      EXPECT_THAT(result.first.restrictedMatrix()(row, col), DoubleNear(J(row, col), 1e-8));
      EXPECT_THAT(result.second.restrictedMatrix()(row, col), DoubleNear(K(row, col), 1e-8));
    }
  }
}

TEST_F(FockMatrixTest, TestRestrictedJK_CauchySchwartz_unrestricted) {
  omp_set_num_threads(1);

  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto scineAtoms = Utils::XyzStreamHandler::read(xyzInput);

  // Reference
  Eigen::MatrixXd D_alpha;
  D_alpha.resize(11, 11);
  Eigen::MatrixXd D_beta;
  D_beta.resize(11, 11);
  Eigen::MatrixXd J_alpha;
  J_alpha.resize(11, 11);
  Eigen::MatrixXd K_alpha;
  K_alpha.resize(11, 11);

  D_alpha << 0.0413176955, -0.0254942060, 0.2171627392, -0.0000000000, 0.0000000000, 0.2741156191, 0.0113652297,
      -0.1945506398, 0.0000000000, 0.0000000000, 0.1204465060, -0.0254942060, 1.0364725940, -0.1892548225, -0.0000000000,
      -0.0000000000, -0.2154639449, -0.0073382827, 0.1502936129, 0.0000000000, 0.0000000000, -0.0758205752,
      0.2171627392, -0.1892548225, 1.4160378550, -0.0000000000, 0.0000000000, 1.6258086504, 0.1471573854, -1.5228384945,
      -0.0000000000, -0.0000000000, 1.1992578529, -0.0000000000, -0.0000000000, -0.0000000000, 0.0001655408,
      0.0000000000, -0.0000000000, 0.0000000000, -0.0000000000, -0.0129115409, -0.0000000000, -0.0000000000,
      0.0000000000, -0.0000000000, 0.0000000000, 0.0000000000, 0.0001655408, 0.0000000000, 0.0000000000, -0.0000000000,
      -0.0000000000, -0.0129115409, 0.0000000000, 0.2741156191, -0.2154639449, 1.6258086504, -0.0000000000,
      0.0000000000, 1.9434902837, 0.1317599770, -1.6182817077, 0.0000000000, 0.0000000000, 1.1857443866, 0.0113652297,
      -0.0073382827, 0.1471573854, 0.0000000000, 0.0000000000, 0.1317599770, 1.0707739569, -0.3806471206, -0.0000000000,
      -0.0000000000, 0.1021855214, -0.1945506398, 0.1502936129, -1.5228384945, -0.0000000000, -0.0000000000,
      -1.6182817077, -0.3806471206, 2.5360360175, 0.0000000000, 0.0000000000, -1.1238882349, 0.0000000000, 0.0000000000,
      -0.0000000000, -0.0129115409, -0.0000000000, 0.0000000000, -0.0000000000, 0.0000000000, 1.0070500047,
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000, -0.0129115409, 0.0000000000,
      -0.0000000000, 0.0000000000, 0.0000000000, 1.0070500047, -0.0000000000, 0.1204465060, -0.0758205752, 1.1992578529,
      -0.0000000000, 0.0000000000, 1.1857443866, 0.1021855214, -1.1238882349, 0.0000000000, -0.0000000000, 1.8616836450;

  D_beta = D_alpha;

  J_alpha << 2.0714309845, 0.3752206719, 1.4213685500, -0.0000000000, 0.0000000000, -1.0579559307, 0.0284170292,
      0.2634298035, -0.0000000000, 0.0000000000, -0.3120957808, 0.3752206719, 6.7358689009, 1.3551494919, -0.0000000000,
      0.0000000000, 0.0893999524, 0.0000554358, 0.3102555857, -0.0000000000, 0.0000000000, -0.5144323275, 1.4213685500,
      1.3551494919, 3.6168273439, -0.0000000000, 0.0000000000, 0.7375522875, 0.3886904545, 1.9908617235, -0.0000000000,
      0.0000000000, -1.8076267950, -0.0000000000, -0.0000000000, -0.0000000000, 3.4812127177, -0.0000000000,
      -0.0000000000, 0.0000000000, -0.0000000000, 1.1131152529, -0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000000, -0.0000000000, 3.4812127177, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000, 1.1131152529,
      -0.0000000000, -1.0579559307, 0.0893999524, 0.7375522875, -0.0000000000, 0.0000000000, 3.8785493863, 0.6571409210,
      2.3386069579, 0.0000000000, 0.0000000000, -1.7088934946, 0.0284170292, 0.0000554358, 0.3886904545, 0.0000000000,
      0.0000000000, 0.6571409210, 9.2032272172, 1.7890651293, 0.0000000000, 0.0000000000, -0.0288900213, 0.2634298035,
      0.3102555857, 1.9908617235, -0.0000000000, 0.0000000000, 2.3386069579, 1.7890651293, 4.8100024915, 0.0000000000,
      0.0000000000, -0.2710327774, -0.0000000000, -0.0000000000, -0.0000000000, 1.1131152529, -0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 4.7539921072, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000000, -0.0000000000, 1.1131152529, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 4.7539921072,
      -0.0000000000, -0.3120957808, -0.5144323275, -1.8076267950, 0.0000000000, -0.0000000000, -1.7088934946,
      -0.0288900213, -0.2710327774, 0.0000000000, -0.0000000000, 4.9105670878;

  K_alpha << 0.0255677106, 0.1724437643, 0.0685929146, -0.0000000000, -0.0000000000, 0.0111331645, 0.0125793071,
      0.0459287305, -0.0000000000, -0.0000000000, -0.0427669367, 0.1724437643, 3.5987745447, 0.6018156001,
      -0.0000000000, -0.0000000000, 0.0134623381, -0.0001366383, 0.1372239892, -0.0000000000, 0.0000000000,
      -0.2277444116, 0.0685929146, 0.6018156001, 0.4896952498, -0.0000000000, 0.0000000000, 0.4537547557, 0.1662559570,
      0.4746985245, -0.0000000000, -0.0000000000, -0.3645407374, -0.0000000000, -0.0000000000, -0.0000000000,
      0.1569121437, -0.0000000000, -0.0000000000, 0.0000000000, -0.0000000000, 0.2167774145, -0.0000000000,
      -0.0000000000, -0.0000000000, -0.0000000000, 0.0000000000, -0.0000000000, 0.1569121437, 0.0000000000,
      -0.0000000000, 0.0000000000, -0.0000000000, 0.2167774145, -0.0000000000, 0.0111331645, 0.0134623381, 0.4537547557,
      -0.0000000000, 0.0000000000, 0.6745130143, 0.2781079508, 0.6350847233, 0.0000000000, 0.0000000000, -0.4212469326,
      0.0125793071, -0.0001366383, 0.1662559570, 0.0000000000, -0.0000000000, 0.2781079508, 4.3010070452, 0.7480005301,
      -0.0000000000, -0.0000000000, 0.0015782095, 0.0459287305, 0.1372239892, 0.4746985245, -0.0000000000, 0.0000000000,
      0.6350847233, 0.7480005301, 1.2832094433, 0.0000000000, -0.0000000000, -0.0475121851, -0.0000000000, -0.0000000000,
      -0.0000000000, 0.2167774145, -0.0000000000, 0.0000000000, -0.0000000000, 0.0000000000, 1.0428330658, 0.0000000000,
      0.0000000000, -0.0000000000, 0.0000000000, -0.0000000000, -0.0000000000, 0.2167774145, 0.0000000000, -0.0000000000,
      -0.0000000000, 0.0000000000, 1.0428330658, 0.0000000000, -0.0427669367, -0.2277444116, -0.3645407374,
      -0.0000000000, -0.0000000000, -0.4212469326, 0.0015782095, -0.0475121851, 0.0000000000, 0.0000000000, 1.1033966945;

  Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(11, 11);
  Eigen::MatrixXd coeffs2 = coeffs;

  Utils::MolecularOrbitals mos = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients(coeffs, coeffs2);

  Utils::LcaoUtils::DensityMatrixBuilder builder(mos);
  auto densityMatrix = builder.generateUnrestrictedForNumberAlphaAndBetaElectrons(7, 7);

  densityMatrix.addMatrixRestricted(D_alpha + D_beta);
  densityMatrix.addMatrixAlpha(D_alpha);
  densityMatrix.addMatrixBeta(D_beta);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", true);

  auto basis = eval.initializeBasisSet("sto-3g", scineAtoms);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;

  auto prescreener = TwoBody::CauchySchwarzDensityPrescreener(basis, densityMatrix);

  std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix> result;

  auto saver = TwoBody::CoulombExchangeDigester(basis, basis, specifier, densityMatrix);
  auto evaluator = TwoBody::Evaluator<TwoBody::CoulombExchangeDigester, TwoBody::CauchySchwarzDensityPrescreener>(
      basis, basis, specifier, std::move(saver), std::move(prescreener));
  evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
  result = evaluator.getResult();

  for (int row = 0; row < 10; ++row) {
    for (int col = 0; col < 10; ++col) {
      EXPECT_THAT(result.first.alphaMatrix()(row, col), DoubleNear(J_alpha(row, col), 1e-8));
      EXPECT_THAT(result.second.alphaMatrix()(row, col), DoubleNear(K_alpha(row, col), 1e-8));
      EXPECT_THAT(result.first.betaMatrix()(row, col), DoubleNear(J_alpha(row, col), 1e-8));
      EXPECT_THAT(result.second.betaMatrix()(row, col), DoubleNear(K_alpha(row, col), 1e-8));
    }
  }
}

TEST_F(FockMatrixTest, TestJK_nuclearElectronic) {
  omp_set_num_threads(1);

  Eigen::MatrixXd L1;
  L1.resize(11, 11);
  Eigen::MatrixXd D1;
  D1.resize(11, 11);
  Eigen::MatrixXd D2;
  D2.resize(5, 5);

  L1 << -0.6988454884, -0.0279938907, -0.2693213832, 0.0000000000, 0.0000000000, 0.3100785137, -0.0008131871,
      -0.0222333572, 0.0000000000, 0.0000000000, 0.0321415240, -0.0279938907, -0.3819127779, -0.0946485662,
      0.0000000000, 0.0000000000, 0.0085872881, -0.0000024761, -0.0199457480, 0.0000000000, 0.0000000000, 0.0330124007,
      -0.2693213832, -0.0946485662, -0.3699706261, 0.0000000000, 0.0000000000, 0.0924686089, -0.0110502969,
      -0.1196704713, 0.0000000000, 0.0000000000, 0.1268569569, 0.0000000000, 0.0000000000, 0.0000000000, -0.3541776036,
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.0731528677, 0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, -0.3541776036, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      -0.0731528677, 0.0000000000, 0.3100785137, 0.0085872881, 0.0924686089, 0.0000000000, 0.0000000000, -0.4009286258,
      -0.0186538488, -0.1100465753, 0.0000000000, 0.0000000000, 0.0827364612, -0.0008131871, -0.0000024761,
      -0.0110502969, 0.0000000000, 0.0000000000, -0.0186538488, -0.2145326945, -0.0504185237, 0.0000000000,
      0.0000000000, 0.0026425285, -0.0222333572, -0.0199457480, -0.1196704713, 0.0000000000, 0.0000000000,
      -0.1100465753, -0.0504185237, -0.2141084396, 0.0000000000, 0.0000000000, 0.0332491721, 0.0000000000, 0.0000000000,
      0.0000000000, -0.0731528677, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.2105256414, 0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.0731528677, 0.0000000000, 0.0000000000,
      0.0000000000, 0.0000000000, -0.2105256414, 0.0000000000, 0.0321415240, 0.0330124007, 0.1268569569, 0.0000000000,
      0.0000000000, 0.0827364612, 0.0026425285, 0.0332491721, 0.0000000000, 0.0000000000, -0.2212480811;

  D1 << 0.0586801148, -0.0461163619, 0.3909275847, 0.0000000000, 0.0000000000, 0.4270165993, 0.0202804578, -0.3249541995,
      0.0000000000, -0.0000000000, 0.1790789829, -0.0461163619, 2.0887351917, -0.4684039750, -0.0000000000,
      -0.0000000000, -0.4571264669, -0.0172661514, 0.3404015775, 0.0000000000, 0.0000000000, -0.1936681390, 0.3909275847,
      -0.4684039750, 3.3182352520, 0.0000000000, 0.0000000000, 3.2931027798, 0.3061307803, -3.2150553786, -0.0000000000,
      -0.0000000000, 2.5809174950, 0.0000000000, -0.0000000000, 0.0000000000, 0.0020052512, -0.0000000000, 0.0000000000,
      -0.0000000000, -0.0000000000, -0.0638595639, 0.0000000000, 0.0000000000, 0.0000000000, -0.0000000000,
      0.0000000000, -0.0000000000, 0.0020052512, 0.0000000000, -0.0000000000, -0.0000000000, 0.0000000000,
      -0.0638595639, 0.0000000000, 0.4270165993, -0.4571264669, 3.2931027798, 0.0000000000, 0.0000000000, 3.3891644407,
      0.2548897285, -3.0262023947, -0.0000000000, -0.0000000000, 2.1703827296, 0.0202804578, -0.0172661514,
      0.3061307803, -0.0000000000, -0.0000000000, 0.2548897285, 2.1416577286, -0.7611365036, -0.0000000000, 0.0000000000,
      0.2046889106, -0.3249541995, 0.3404015775, -3.2150553786, -0.0000000000, -0.0000000000, -3.0262023947,
      -0.7611365036, 5.0305350709, -0.0000000000, 0.0000000000, -2.2144648021, 0.0000000000, 0.0000000000, -0.0000000000,
      -0.0638595639, 0.0000000000, -0.0000000000, -0.0000000000, -0.0000000000, 2.0336823242, 0.0000000000,
      -0.0000000000, -0.0000000000, 0.0000000000, -0.0000000000, 0.0000000000, -0.0638595639, -0.0000000000,
      0.0000000000, 0.0000000000, 0.0000000000, 2.0336823242, -0.0000000000, 0.1790789829, -0.1936681390, 2.5809174950,
      0.0000000000, 0.0000000000, 2.1703827296, 0.2046889106, -2.2144648021, -0.0000000000, -0.0000000000, 3.6984800533;

  D2 << 0.0497317488, 0.1438798072, -0.0000000000, -0.0000000000, -0.1309061669, 0.1438798072, 0.4162612297,
      -0.0000000000, 0.0000000000, -0.3787269602, 0.0000000000, -0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000000, -0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, -0.1309061669, -0.3787269602, 0.0000000000,
      0.0000000000, 0.3445771553;

  Eigen::MatrixXd coeffs = Eigen::MatrixXd::Random(11, 11);
  Utils::MolecularOrbitals mos = Utils::MolecularOrbitals::createFromRestrictedCoefficients(coeffs);
  Utils::LcaoUtils::DensityMatrixBuilder builder(mos);
  auto densityMatrix1 = builder.generateRestrictedForNumberElectrons(14);
  // densityMatrix1.restrictedMatrix() = D1;
  Eigen::MatrixXd D1cp = D1;
  densityMatrix1.setDensity(std::move(D1cp), 14);

  Eigen::MatrixXd coeffs1 = Eigen::MatrixXd::Random(5, 5);
  Eigen::MatrixXd coeffs2 = Eigen::MatrixXd::Random(5, 5);
  Utils::MolecularOrbitals mos_p = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients(coeffs1, coeffs2);
  Utils::LcaoUtils::DensityMatrixBuilder builder_p(mos_p);
  Utils::DensityMatrix densityMatrix2 = builder_p.generateUnrestrictedForNumberAlphaAndBetaElectrons(1, 0);
  Eigen::MatrixXd D2_alpha = D2;
  Eigen::MatrixXd D2_beta;
  D2_beta.setZero(D2.rows(), D2.cols());
  densityMatrix2.setDensity(std::move(D2_alpha), std::move(D2_beta), 1, 0);

  std::stringstream xyzInput_e("3\n\n"
                               "H    0.0     0.0 0.0\n"
                               "C    1.0378  0.0 0.0\n"
                               "N    2.204   0.0 0.0");

  auto atoms1 = Utils::XyzStreamHandler::read(xyzInput_e);

  std::stringstream xyzInput_p("1\n\n"
                               "H    0.0     0.0 0.0");

  auto atoms2 = Utils::XyzStreamHandler::read(xyzInput_p);

  LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", true);

  auto basis1 = eval.initializeBasisSet("sto-3g", atoms1);
  auto basis2 = eval.initializeBasisSet("6-31g**", atoms2);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;
  specifier.typeVector = {Utils::Integrals::getParticleType("e"), Utils::Integrals::getParticleType("H")};

  auto prescreener = Integrals::TwoBody::VoidPrescreener();

  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> result;

  auto saver_direct = Integrals::TwoBody::TwoTypeCoulombDigester(basis1, basis2, specifier, densityMatrix1, densityMatrix2);
  auto evaluator_direct = Integrals::TwoBody::Evaluator<Integrals::TwoBody::TwoTypeCoulombDigester>(
      basis1, basis2, specifier, std::move(saver_direct), std::move(prescreener));
  evaluator_direct.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
  result = evaluator_direct.getResult();

  auto result_map = LibintIntegrals::evaluate(specifier, basis1, basis2);

  const auto& result_conv = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  auto test1 = result.first;

  std::cout << std::fixed << std::setprecision(10);
  // Ref:
  {
    Eigen::MatrixXd ref(D1.rows(), D1.rows());
    ref.setZero();
    int d1 = int(D1.rows());
    int d2 = int(D2.rows());

    for (int i = 0; i < d1; ++i) {
      for (int j = 0; j < d1; ++j) {
        for (int K = 0; K < d2; ++K) {
          for (int L = 0; L < d2; ++L) {
            auto index = Scine::Integrals::TwoBody::IndexType<4>({i, j, K, L});

            double integralValue = result_conv(i * d1 + j, K * d2 + L);

            if (!std::equal(index.begin(), index.end(),
                            Integrals::TwoBody::getMappedIndex<Integrals::TwoBody::IntegralSymmetry::fourfold>(index).begin())) {
              continue;
            }

            auto symmetricIndexes =
                Scine::Integrals::TwoBody::getSymmetricIndices<Scine::Integrals::TwoBody::IntegralSymmetry::fourfold>(
                    std::move(index));
            // Coulomb term
            for (const auto& idx : symmetricIndexes) {
              ref(idx[0], idx[1]) += integralValue * D2(idx[2], idx[3]);
            }
          }
        }
      }
    }

    for (int row = 0; row < test1.rows(); ++row) {
      for (int col = 0; col < test1.cols(); ++col) {
        EXPECT_THAT(ref(row, col), DoubleNear(L1(row, col), 1e-8));
      }
    }
  }

  for (int row = 0; row < test1.rows(); ++row) {
    for (int col = 0; col < test1.cols(); ++col) {
      EXPECT_THAT(test1(row, col), DoubleNear(L1(row, col), 1e-8));
    }
  }
}

// TEST_F(FockMatrixTest, MakeRefernceData) {
//  // Reference
//
//  std::stringstream xyzInput("3\n\n"
//                             "O          0.00000       -0.07579        0.00000\n"
//                             "H          0.86681        0.60144        0.00000\n"
//                             "H         -0.86681        0.60144        0.00000");
//
//  auto scineAtoms = Utils::XyzStreamHandler::read(xyzInput);
//
//  LibintIntegrals eval;
//
//  Utils::Settings& settings = eval.settings();
//
//  settings.modifyBool("use_pure_spherical", true);
//
//  auto basis = eval.initializeBasisSet("def2-svp", scineAtoms);
//
//  Utils::Integrals::IntegralSpecifier specifier;
//  specifier.op = Utils::Integrals::Operator::Overlap;
//
//  auto result_map_S = LibintIntegrals::evaluate(specifier, basis, basis);
//
//  const auto& result_S = result_map_S[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
//
//  Utils::Integrals::IntegralSpecifier point_charges;
//  point_charges.op = Utils::Integrals::Operator::PointCharges;
//  point_charges.atoms = scineAtoms;
//
//  auto result_map_Vcore = LibintIntegrals::evaluate(point_charges, basis, basis);
//  const auto& result_Vcore = result_map_Vcore[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value,
//  0}];
//
//  Utils::Integrals::IntegralSpecifier kinetic;
//  kinetic.op = Utils::Integrals::Operator::Kinetic;
//
//  auto result_map_T = LibintIntegrals::evaluate(kinetic, basis, basis);
//  const auto& result_T = result_map_T[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
//
//  Eigen::MatrixXd Hcore = result_T + result_Vcore;
//
//  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Hcore, result_S);
//
//  auto coeffs = ges.eigenvectors();
//  Utils::MolecularOrbitals mos = Utils::MolecularOrbitals::createFromRestrictedCoefficients(coeffs);
//  Utils::LcaoUtils::DensityMatrixBuilder builder(mos);
//  auto densityMatrix = builder.generateRestrictedForNumberElectrons(10);
//
//  // Fock matrix
//
//  Utils::Integrals::IntegralSpecifier specifier_jk;
//  specifier_jk.op = Utils::Integrals::Operator::Coulomb;
//
//  auto prescreener = TwoBody::CauchySchwarzDensityPrescreener(basis, densityMatrix);
//
//  std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix> result;
//
//  auto saver = TwoBody::CoulombExchangeDigester(basis, basis, specifier_jk, densityMatrix);
//  auto evaluator = TwoBody::Evaluator<TwoBody::CoulombExchangeDigester, TwoBody::CauchySchwarzDensityPrescreener>(
//      basis, basis, specifier, std::move(saver), std::move(prescreener));
//  evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
//  result = evaluator.getResult();
//
//  Eigen::MatrixXd F = result_T + result_Vcore + result.first.restrictedMatrix() - 0.5 *
//  result.second.restrictedMatrix();
//
//  auto print = [&](const Eigen::MatrixXd mat) {
//    std::cout << std::setprecision(10) << std::fixed;
//    for (auto i = 0L; i < mat.rows(); ++i) {
//      for (auto j = 0L; j < mat.rows(); ++j) {
//        std::cout << std::setw(16) << mat(i, j) << ",";
//      }
//    }
//    std::cout << "\n";
//  };
//
//  std::cout << "S << ";
//  print(result_S);
//  std::cout << "T << ";
//  print(result_T);
//  std::cout << "V << ";
//  print(result_Vcore);
//  std::cout << "F << ";
//  print(F);
//  std::cout << "D << ";
//  print(densityMatrix.restrictedMatrix());
//}
