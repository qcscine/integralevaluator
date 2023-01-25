/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INTEGRALEVALUATOR_DIGESTER_H
#define INTEGRALEVALUATOR_DIGESTER_H

#include <LibintIntegrals/TwoBodyIntegrals/SymmetryHelper.h>
#include <Utils/DataStructures/IntegralSpecifier.h>

namespace Scine {
namespace Utils {
namespace Integrals {
class BasisSet;
}
} // namespace Utils

namespace Integrals {
namespace TwoBody {
/**
 * @class TwoBodiesIntegralsDigester @file TwoBodiesIntegralsDigester
 * @brief CRTP functor representing a static interface for a digester.
 * @tparam CrtpClass The CRTP parameter, i.e. the derived class implementing the interface.
 */
template<typename CrtpClass>
class Digester {
 public:
  using ResultType = void;
  Digester(const Utils::Integrals::BasisSet& scineBasis, const Utils::Integrals::IntegralSpecifier& specifier);
  Digester(const Utils::Integrals::BasisSet& scineBasis1, const Utils::Integrals::BasisSet& scineBasis2,
           const Utils::Integrals::IntegralSpecifier& specifier);
  /**
   * @brief Function to be called in the usage of a digester.
   * Every digester needs to unpack the integral over shells in integral over basis functions. This is why this id done
   * by the TwoBodiesIntegralsDigester. What this function does:
   * - It takes a mapped quartet integral shell as arguments. Depending on the derivative type, it looks up how many
   * different integral types there are.
   * - It loops over the different integral derivative types (derivative of the first shell wrt its nuclear geometric
   * derivative,...) and unpacks the derivative types into the integral over basis functions.
   * - It forwards the integral over basis functions
   *
   * NB: The integrals are scaled by the degeneracy factor. Normally this implies eight-fold symmetry. In order to
   * bypass this default behaviour (i.e. AO2MO), hide the calculateDegeneracy() symbol in the derived class.
   * @param buffer Eigen::MatrixXd containing the shell-quartet integrals. Each row is a different type of integral,
   * i.e. row1 -> first derivative of the first shell with respect to a nuclear geometric coordinate, row2 -> first
   * derivative of the second shell with respect to a nuclear geometric coordinate,...
   */
  void operator()(const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>& buffer, int shell1, int shell2, int shell3, int shell4) {
    double degeneracy = derived().computeDegeneracyImpl(shell1, shell2, shell3, shell4);

    auto performLoop = [&](const int center = 0, const Utils::Integrals::DerivKey& key = Utils::Integrals::DerivKey::value) {
      // Retrieve the correct result by combining center and deriv key
      auto index = center * numRelevantDerivKeys_ + static_cast<int>(key);
      const Eigen::VectorXd& vectorBuffer = buffer.row(index);

      size_t bufferIndex = 0;

      for (size_t functionInShell1 = 0; functionInShell1 < scineBasis1_[shell1].size(); ++functionInShell1) {
        const size_t basisFunction1 = functionInShell1 + indexFirstBFInShell1_[shell1];
        for (size_t functionInShell2 = 0; functionInShell2 < scineBasis1_[shell2].size(); ++functionInShell2) {
          const size_t basisFunction2 = functionInShell2 + indexFirstBFInShell1_[shell2];
          for (size_t functionInShell3 = 0; functionInShell3 < scineBasis2_[shell3].size(); ++functionInShell3) {
            const size_t basisFunction3 = functionInShell3 + indexFirstBFInShell2_[shell3];
            for (size_t functionInShell4 = 0; functionInShell4 < scineBasis2_[shell4].size(); ++functionInShell4) {
              const size_t basisFunction4 = functionInShell4 + indexFirstBFInShell2_[shell4];

              auto const integral = vectorBuffer(bufferIndex);
              ++bufferIndex;
              if (integral != 0.0) {
                derived().digestImpl(integral, basisFunction1, basisFunction2, basisFunction3, basisFunction4, index, degeneracy);
              }
            }
          }
        }
      }
    };
    if (specifier_.derivOrder == 0) {
      performLoop();
    }
    // If derivative, calculate derivative
    else if (specifier_.derivOrder == 1) {
      for (int center = 0; center < numberOfCenters_; ++center) {
        for (auto const& derivKey :
             {Utils::Integrals::DerivKey::x, Utils::Integrals::DerivKey::y, Utils::Integrals::DerivKey::z}) {
          performLoop(center, derivKey);
        }
      }
    }
    assert(specifier_.derivOrder < 2);
  }

  /**
   * @brief Interface method to compute the degeneracy number of the shell quartet. The integrals need to be scaled by
   * the degeneracy. This method computes it. It allows to perform general scaling of the shell integrals. If the
   * underlying digester does not hide the computeDegenerayImpl(int, int, int, int) method, the degeneracy is defaulted
   * to the one given by the 1-fold symmetric eri tensor, so to 1.0. See getDegeneracy<IntegralSymmetry
   * symmetry>(shell1, shell2, shell3, shell4) for details.
   *
   * @param shell1
   * @param shell2
   * @param shell3
   * @param shell4
   * @return The degeneracy of the shell quartet.
   */
  double computeDegeneracy(int shell1, int shell2, int shell3, int shell4) {
    return derived().computeDegeneracyImpl(shell1, shell2, shell3, shell4);
  }
  /**
   * @brief Generic result return function. The type depends on the kind of digester. This is a >= C++14 specific
   * construct. In CRTP, the template parameter is incomplete at this point of compilation. This means that a construct
   * of the type
   * @code{cpp}
   * const typedef CrtpClass::ResultType& getResult() const {
   * // ...
   * }
   * @endcode
   * results in an incomplete-type compilation error. This auto construct delays the type deduction until CrtpClass is
   * complete.
   */
  const auto& getResult() const {
    return derived().getResultImpl();
  }

  /**
   * @brief initializer for the digester. Needed for instance to initialize local matrices for parallelization.
   * @param numberThreads The number of threads with which to parallelize.
   */
  void initialize(int numberThreads) {
    derived().initializeImpl(numberThreads);
  }

  /**
   * @brief Class to allow for operation that have to be done just one time before the result is given back. This class
   * exist in order not to mess with the const-correctness of getResult(). It allows, for example, to symmetrize the
   * 2-electron matrix just one time at the very end in such a way that it is no more needed to symmetrize every single
   * operation in the digester.
   */
  void finalize() {
    derived().finalizeImpl();
  }

 private:
  inline double computeDegeneracyImpl(int shell1, int shell2, int shell3, int shell4) const {
    return getDegeneracy<IntegralSymmetry::onefold>(shell1, shell2, shell3, shell4);
  }

  CrtpClass& derived() {
    return static_cast<CrtpClass&>(*this);
  }

  const CrtpClass& derived() const {
    return static_cast<const CrtpClass&>(*this);
  }

 protected:
  const Utils::Integrals::BasisSet& scineBasis1_;
  const Utils::Integrals::BasisSet& scineBasis2_;
  const Utils::Integrals::IntegralSpecifier& specifier_;
  const int numberOfCenters_;
  const std::vector<size_t> indexFirstBFInShell1_;
  const std::vector<size_t> indexFirstBFInShell2_;
  const std::size_t dim1_;
  const std::size_t dim2_;
  const std::size_t dim1sq_;
  const std::size_t dim2sq_;
  const int numRelevantDerivKeys_ = 3;
  double scaling_ = 1.;
};

template<typename CrtpClass>
Digester<CrtpClass>::Digester(const Utils::Integrals::BasisSet& scineBasis1, const Utils::Integrals::BasisSet& scineBasis2,
                              const Utils::Integrals::IntegralSpecifier& specifier)
  : scineBasis1_(scineBasis1),
    scineBasis2_(scineBasis2),
    specifier_(specifier),
    numberOfCenters_((specifier_.derivOrder > 0) ? 4 : 1),
    indexFirstBFInShell1_(scineBasis1_.shell2bf()),
    indexFirstBFInShell2_(scineBasis2_.shell2bf()),
    dim1_(scineBasis1.nbf()),
    dim2_(scineBasis2.nbf()),
    dim1sq_(dim1_ * dim1_),
    dim2sq_(dim2_ * dim2_) {
}

template<typename CrtpClass>
Digester<CrtpClass>::Digester(const Utils::Integrals::BasisSet& scineBasis, const Utils::Integrals::IntegralSpecifier& specifier) {
  Digester(scineBasis, scineBasis, specifier);
}

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_DIGESTER_H
