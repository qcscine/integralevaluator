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

using namespace Scine;
using namespace Integrals;

using namespace testing;

class BasisSetTest : public Test {};

class BasisSetHandlerTest : public Test {};

class ShellPairTest : public Test {};

TEST_F(BasisSetTest, shell2atom) {
  std::stringstream Ethanol("9\n\n"
                            "C    -4.0410150   -1.2118929   -0.0394793 \n"
                            "C    -2.6937552   -0.5748019    0.1881221 \n"
                            "O    -1.7026307   -1.5718950    0.0262698 \n"
                            "H    -4.8396201   -0.4780258    0.0891541 \n"
                            "H    -4.2018174   -2.0266305    0.6693109 \n"
                            "H    -4.1063106   -1.6202046   -1.0499636 \n"
                            "H    -2.5412751    0.2452513   -0.5277968 \n"
                            "H    -2.6490005   -0.1433850    1.1978935 \n"
                            "H    -0.8426606   -1.1792649    0.1904767 ");

  auto scineAtoms = Utils::XyzStreamHandler::read(Ethanol);
  auto libintAtoms = BasisSetHandler::scineToLibint(scineAtoms);

  LibintIntegrals eval;

  auto scineBasis = eval.initializeBasisSet("cc-pvqz", scineAtoms);

  auto libintBasis = libint2::BasisSet("cc-pvqz", libintAtoms);

  auto libintS2A = libintBasis.shell2atom(libintAtoms);

  auto scineS2A = scineBasis.shellToAtom(scineAtoms);

  for (std::size_t i = 0; i < libintS2A.size(); ++i) {
    ASSERT_EQ(libintS2A.at(i), scineS2A.at(i));
  }
}

TEST_F(BasisSetTest, atom2shell) {
  std::stringstream Ethanol("9\n\n"
                            "C    -4.0410150   -1.2118929   -0.0394793 \n"
                            "C    -2.6937552   -0.5748019    0.1881221 \n"
                            "O    -1.7026307   -1.5718950    0.0262698 \n"
                            "H    -4.8396201   -0.4780258    0.0891541 \n"
                            "H    -4.2018174   -2.0266305    0.6693109 \n"
                            "H    -4.1063106   -1.6202046   -1.0499636 \n"
                            "H    -2.5412751    0.2452513   -0.5277968 \n"
                            "H    -2.6490005   -0.1433850    1.1978935 \n"
                            "H    -0.8426606   -1.1792649    0.1904767 ");

  auto scineAtoms = Utils::XyzStreamHandler::read(Ethanol);
  auto libintAtoms = BasisSetHandler::scineToLibint(scineAtoms);

  LibintIntegrals eval;

  auto scineBasis = eval.initializeBasisSet("cc-pvqz", scineAtoms);

  auto libintBasis = libint2::BasisSet("cc-pvqz", libintAtoms);

  auto libintA2S = libintBasis.atom2shell(libintAtoms);

  auto scineA2S = scineBasis.atomToShell(scineAtoms);

  for (std::size_t i = 0; i < libintA2S.size(); ++i) {
    for (std::size_t j = 0; j < libintA2S[i].size(); ++j) {
      ASSERT_EQ(libintA2S[i][j], scineA2S[i][j]);
    }
  }
}

TEST_F(BasisSetTest, DifferentBasisSets) {
  std::stringstream Ethanol("9\n\n"
                            "H    -4.8396201   -0.4780258    0.0891541 \n"
                            "H    -4.2018174   -2.0266305    0.6693109 \n"
                            "H    -4.1063106   -1.6202046   -1.0499636 \n"
                            "H    -2.5412751    0.2452513   -0.5277968 \n"
                            "H    -2.6490005   -0.1433850    1.1978935 \n"
                            "H    -0.8426606   -1.1792649    0.1904767 \n"
                            "C    -4.0410150   -1.2118929   -0.0394793 \n"
                            "C    -2.6937552   -0.5748019    0.1881221 \n"
                            "O    -1.7026307   -1.5718950    0.0262698 \n");

  auto scineAtoms = Utils::XyzStreamHandler::read(Ethanol);
  auto libintAtoms = BasisSetHandler::scineToLibint(scineAtoms);

  LibintIntegrals eval;

  auto basis = eval.initializeBasisSet(
      {{Utils::ElementType::H, "def2-tzvp"}, {Utils::ElementType::C, "aug-cc-pvtz"}, {Utils::ElementType::O, "cc-pvqz"}},
      scineAtoms);

  auto count = 0;
  for (auto i = 0; i < scineAtoms.size(); ++i) {
    for (auto j = 0; j < basis.getAtoms().size(); ++j) {
      if (scineAtoms[i].getPosition() == basis.getAtoms()[j].getPosition())
        ++count;
    }
  }

  ASSERT_TRUE(count == 9);
}

TEST_F(ShellPairTest, CanAddAShellPair) {
  std::stringstream h2("2\n\n"
                       "H 0 0 0\n"
                       "H 1.2 0 0");
  auto scineAtoms = Utils::XyzStreamHandler::read(h2);

  LibintIntegrals eval;

  auto scineBasis = eval.initializeBasisSet("def2-svp", scineAtoms);

  auto libintShell1 = BasisSetHandler::scineToLibint(scineBasis[0]);
  auto libintShell2 = BasisSetHandler::scineToLibint(scineBasis[1]);
  Utils::Integrals::ShellPair shellPair(scineBasis[0]);
  BasisSetHandler::addPair(shellPair, 1, scineBasis[1]);

  Utils::Integrals::ShellPairs shellPairs;
  shellPairs.emplace_back(std::move(shellPair));
  ASSERT_EQ(shellPairs[0].size(), 1);
}

TEST_F(BasisSetHandlerTest, AtomsScineToLibint) {
  Utils::ElementTypeCollection etc{Utils::ElementType::H, Utils::ElementType::H, Utils::ElementType::O};
  Utils::PositionCollection pc = Eigen::MatrixX3d::Random(3, 3);
  Utils::AtomCollection scineAtoms(etc, pc);

  std::vector<libint2::Atom> libintAtoms = BasisSetHandler::scineToLibint(scineAtoms);
  for (auto i = 0; i < scineAtoms.size(); ++i) {
    ASSERT_EQ(static_cast<int>(Utils::ElementInfo::Z(scineAtoms.getElement(i))), libintAtoms.at(i).atomic_number);
    ASSERT_EQ(scineAtoms.getPosition(i)(0), libintAtoms.at(i).x);
    ASSERT_EQ(scineAtoms.getPosition(i)(1), libintAtoms.at(i).y);
    ASSERT_EQ(scineAtoms.getPosition(i)(2), libintAtoms.at(i).z);
  }
}

TEST_F(BasisSetHandlerTest, TestBasisSet) {
  Utils::ElementTypeCollection etc{Utils::ElementType::H, Utils::ElementType::H, Utils::ElementType::O};
  Utils::PositionCollection pc = Eigen::MatrixX3d::Random(3, 3);
  Utils::AtomCollection scineAtoms(etc, pc);

  LibintIntegrals eval;

  auto scineBasis = eval.initializeBasisSet("6-31G*", scineAtoms);

  auto libintAtoms = BasisSetHandler::scineToLibint(scineAtoms);

  auto libintBasis = libint2::BasisSet("6-31G*", libintAtoms);

  libintBasis.set_pure(true);

  ASSERT_EQ(libintBasis.max_l(), scineBasis.max_l());
  ASSERT_EQ(libintBasis.max_nprim(), scineBasis.max_nprim());
  ASSERT_EQ(libintBasis.nbf(), scineBasis.nbf());

  int count = 0;
  for (std::size_t i = 0; i < libintBasis.size(); ++i) {
    // exponents
    for (std::size_t k = 0; k < libintBasis.at(i).alpha.size(); ++k) {
      ASSERT_DOUBLE_EQ(libintBasis.at(i).alpha.at(k), scineBasis.at(i).getVecAlpha().at(k));
    }
    // shift
    for (std::size_t k = 0; k < libintBasis.at(i).O.size(); ++k) {
      ASSERT_DOUBLE_EQ(libintBasis.at(i).O.at(k), scineBasis.at(i).getShift()(k));
    }
    for (std::size_t j = 0; j < libintBasis.at(i).contr.size(); ++j) {
      const auto& contr = libintBasis.at(i).contr.at(j);
      // angular momentum
      ASSERT_EQ(contr.l, scineBasis.at(count).l());
      // contraction
      ASSERT_EQ(contr.coeff.size(), scineBasis.at(count).getContrLegnth());
      for (std::size_t k = 0; k < contr.coeff.size(); ++k) {
        ASSERT_DOUBLE_EQ(contr.coeff.at(k), scineBasis.at(count).getVecCoeffs().at(k));
      }
      ++count;
    }
  }
}
