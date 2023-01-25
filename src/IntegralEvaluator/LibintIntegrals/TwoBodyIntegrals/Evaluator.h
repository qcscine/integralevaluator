/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_EVALUATOR_H
#define INTEGRALEVALUATOR_EVALUATOR_H

#include <LibintIntegrals/BasisSetHandler.h>
#include <LibintIntegrals/TwoBodyIntegrals/VoidPrescreener.h>
#include <memory>

namespace Scine {
namespace Integrals {
namespace TwoBody {

template<typename DigesterType, typename PrescreenerType = VoidPrescreener>
class Evaluator {
 public:
  Evaluator(const Utils::Integrals::BasisSet& scineBasis, const Utils::Integrals::IntegralSpecifier& specifier,
            DigesterType&& digester, PrescreenerType&& prescreener);

  Evaluator(const Utils::Integrals::BasisSet& scineBasis1, const Utils::Integrals::BasisSet& scineBasis2,
            const Utils::Integrals::IntegralSpecifier& specifier, DigesterType&& digester, PrescreenerType&& prescreener);

  const auto& getResult() const {
    return digester_.getResult();
  }

  template<libint2::Operator op>
  void evaluateTwoBodyIntegrals() {
    std::shared_ptr<Utils::Integrals::ShellPairs> shellPairs1 = scineBasis1_.getShellPairs();
    std::shared_ptr<Utils::Integrals::ShellPairs> shellPairs2 = scineBasis2_.getShellPairs();

    Libint::getInstance();
    digester_.initialize(Libint::getMaxNumberThreads());

    // 4 centers times 3 coordinates
    auto const numberOfResults = (specifier_.derivOrder > 0) ? 4 * 3 : 1;

#pragma omp parallel
    {
      auto localEngine = Libint::getEngine(scineBasis1_, scineBasis2_, op, specifier_.derivOrder);
      auto const& buffer = localEngine.results();

      // make libint basis:
      std::vector<libint2::Shell> libintShellVector1;
      libintShellVector1.reserve(scineBasis1_.size());
      for (auto const& shell : scineBasis1_) {
        libintShellVector1.push_back(BasisSetHandler::scineToLibint(shell));
      }
      std::vector<std::vector<libint2::ShellPair>> libintPreComputShellPairVector1;
      libintPreComputShellPairVector1.reserve(shellPairs1->size());
      std::vector<libint2::Shell> libintShellPairVector1;
      libintShellPairVector1.reserve(shellPairs1->size());
      for (auto s = 0; s < shellPairs1->size(); ++s) {
        libintShellPairVector1.push_back(BasisSetHandler::scineToLibint(shellPairs1->at(s).getShell()));
        std::vector<libint2::ShellPair> tmp;
        tmp.reserve(shellPairs1->size());
        for (const Utils::Integrals::ShellPairData& shellPair : shellPairs1->at(s)) {
          tmp.push_back(BasisSetHandler::scineToLibint(*shellPair.precomputedShellPair));
        }
        libintPreComputShellPairVector1.push_back(tmp);
      }
      std::vector<libint2::Shell> libintShellVector2;
      libintShellVector2.reserve(scineBasis2_.size());
      for (auto const& shell : scineBasis2_) {
        libintShellVector2.push_back(BasisSetHandler::scineToLibint(shell));
      }
      std::vector<std::vector<libint2::ShellPair>> libintPreComputShellPairVector2;
      libintPreComputShellPairVector2.reserve(shellPairs2->size());
      std::vector<libint2::Shell> libintShellPairVector2;
      libintShellPairVector2.reserve(shellPairs2->size());
      for (auto s = 0; s < shellPairs2->size(); ++s) {
        libintShellPairVector2.push_back(BasisSetHandler::scineToLibint(shellPairs2->at(s).getShell()));
        std::vector<libint2::ShellPair> tmp;
        tmp.reserve(shellPairs2->size());
        for (const Utils::Integrals::ShellPairData& shellPair : shellPairs2->at(s)) {
          tmp.push_back(BasisSetHandler::scineToLibint(*shellPair.precomputedShellPair));
        }
        libintPreComputShellPairVector2.push_back(tmp);
      }

#pragma omp for schedule(dynamic)
      for (auto s1 = 0UL; s1 < scineBasis1_.size(); ++s1) {
        const auto shell1Size = scineBasis1_[s1].size();
        auto const& pairsOfShell1 = shellPairs1->at(s1);

        for (auto sp12 = 0; sp12 < pairsOfShell1.size(); ++sp12) {
          auto const& shellPair12 = pairsOfShell1[sp12];
          const auto shell2Size = scineBasis1_.at(shellPair12.secondShellIndex).size();
          auto s3_max = (scineBasis1_ == scineBasis2_) ? s1 : scineBasis2_.size() - 1;
          // Account for two-fold symmetry.
          for (auto s3 = 0UL; s3 <= s3_max; ++s3) {
            const auto shell3Size = scineBasis2_[s3].size();
            auto const& pairsOfShell3 = shellPairs2->at(s3);

            int sp34_max;
            if (scineBasis1_ != scineBasis2_) {
              sp34_max = pairsOfShell3.size() - 1;
            }
            else {
              // TODO -> account for 4-fold symmetric of integrals. This assumes that at least to some degree there is
              // an 8-fold symmetry
              //         i.e. even though the COM integrals are 4-fold symmetric, still only the 8-fold symmetric
              //         integrals are required.
              // The reason for this is that the shell pairs are screened. Hence, some can be discarded which
              // must be taken into account.
              if (int(s3) < pairsOfShell3.size()) {
                sp34_max = (s1 == s3) ? sp12 : s3;
              }
              else {
                sp34_max = (s1 == s3) ? sp12 : pairsOfShell3.size() - 1;
              }
            }
            for (auto sp34 = 0; sp34 <= sp34_max; ++sp34) {
              auto const& shellPair34 = pairsOfShell3[sp34];
              const auto shell4Size = scineBasis2_[shellPair34.secondShellIndex].size();
              bool isSignificant = prescreener_(s1, shellPair12.secondShellIndex, s3, shellPair34.secondShellIndex,
                                                shellPair12.cauchySchwarzFactor * shellPair34.cauchySchwarzFactor);

              if (!isSignificant) {
                continue;
              }
              auto& shell1 = libintShellVector1[shellPair12.secondShellIndex];
              const auto* ptrlibintShellPair12 = &libintPreComputShellPairVector1[s1][sp12];
              const auto* ptrlibintShellPair34 = &libintPreComputShellPairVector2[s3][sp34];
              auto& shell3 = libintShellVector2[shellPair34.secondShellIndex];
              if (this->specifier_.derivOrder == 0) {
                localEngine.template compute2<op, libint2::BraKet::xx_xx, static_cast<std::size_t>(0)>(
                    libintShellPairVector1.at(s1), shell1, libintShellPairVector2.at(s3), shell3, ptrlibintShellPair12,
                    ptrlibintShellPair34);
              }
              else if (this->specifier_.derivOrder == 1) {
                localEngine.template compute2<op, libint2::BraKet::xx_xx, static_cast<std::size_t>(1)>(
                    libintShellPairVector1.at(s1), shell1, libintShellPairVector2.at(s3), shell3, ptrlibintShellPair12,
                    ptrlibintShellPair34);
              }
              /* Everything is screened out */
              if (buffer[0] == nullptr) {
                continue;
              }

              Eigen::Matrix<double, -1, -1, Eigen::RowMajor> bufferMatrix(numberOfResults, shell1Size * shell2Size *
                                                                                               shell3Size * shell4Size);
              for (int res = 0; res < numberOfResults; ++res) {
                bufferMatrix.row(res) =
                    Eigen::Map<const Eigen::RowVectorXd>(buffer[res], shell1Size * shell2Size * shell3Size * shell4Size);
              }
              digester_(bufferMatrix, s1, shellPair12.secondShellIndex, s3, shellPair34.secondShellIndex);
            }
          }
        }
      }
    }
    digester_.finalize();
  }

 private:
  const Utils::Integrals::BasisSet& scineBasis1_;
  const Utils::Integrals::BasisSet& scineBasis2_;
  const Utils::Integrals::IntegralSpecifier& specifier_;
  DigesterType digester_;
  PrescreenerType prescreener_;
};

template<typename DigesterType, typename PrescreenerType>
Evaluator<DigesterType, PrescreenerType>::Evaluator(const Utils::Integrals::BasisSet& scineBasis,
                                                    const Utils::Integrals::IntegralSpecifier& specifier,
                                                    DigesterType&& digester, PrescreenerType&& prescreener)
  : scineBasis1_(scineBasis), scineBasis2_(scineBasis), digester_(std::move(digester)), prescreener_(std::move(prescreener)) {
}

template<typename DigesterType, typename PrescreenerType>
Evaluator<DigesterType, PrescreenerType>::Evaluator(const Utils::Integrals::BasisSet& scineBasis1,
                                                    const Utils::Integrals::BasisSet& scineBasis2,
                                                    const Utils::Integrals::IntegralSpecifier& specifier,
                                                    DigesterType&& digester, PrescreenerType&& prescreener)
  : scineBasis1_(scineBasis1),
    scineBasis2_(scineBasis2),
    specifier_(specifier),
    digester_(std::move(digester)),
    prescreener_(std::move(prescreener)) {
}

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_EVALUATOR_H
