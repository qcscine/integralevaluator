/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INTEGRALEVALUATOR_COULOMBEXCHANGECONSTRUCTOR_H
#define INTEGRALEVALUATOR_COULOMBEXCHANGECONSTRUCTOR_H

#include <Utils/DataStructures/SpinAdaptedMatrix.h>

namespace Scine {
namespace Utils {
class DensityMatrix;
}
namespace Integrals {
namespace TwoBody {

/**
 * @class TwoElectronMatrixEvaluator @file TwoElectronMatrixEvaluator.h
 * @brief This class evaluates the two electron matrix elements.
 * This class evaluates the two electron matrix by evaluating integrals of
 * basis functions quartet and contracting them with the density.
 * It avoids having to duplicate the integral to matrix code every instance
 * it comes up.
 */
class CoulombExchangeConstructor {
 public:
  explicit CoulombExchangeConstructor(const Utils::DensityMatrix& densityMatrix);

  /**
   * @brief Evaluates 6 matrix elements from a basis function quartet.
   */
  void evaluateBasisFunctionQuartet(double integralValue, int basisFunction1, int basisFunction2, int basisFunction3,
                                    int basisFunction4);

  /**
   * @brief Finalized the calculation, i.e. symmetrize the result.
   */
  void finalizeEvaluation();

  /**
   * @brief Getter for the Coulomb matrix
   */
  const Utils::SpinAdaptedMatrix& getCoulombMatrix() const;
  /**
   * @brief Getter for the exchange matrix
   */
  const Utils::SpinAdaptedMatrix& getExchangeMatrix() const;

 private:
  Utils::SpinAdaptedMatrix coulomb_;
  Utils::SpinAdaptedMatrix exchange_;
  const Utils::DensityMatrix& densityMatrix_;
  unsigned long dim_;
};

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_COULOMBEXCHANGECONSTRUCTOR_H
