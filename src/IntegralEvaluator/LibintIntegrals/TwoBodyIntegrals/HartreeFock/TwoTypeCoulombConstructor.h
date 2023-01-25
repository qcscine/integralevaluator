/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INTEGRALEVALUATOR_TWOTYPECOULOMBCONSTRUCTOR_H
#define INTEGRALEVALUATOR_TWOTYPECOULOMBCONSTRUCTOR_H

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
class TwoTypeCoulombConstructor {
 public:
  explicit TwoTypeCoulombConstructor(const Utils::DensityMatrix& densityMatrix_type1,
                                     const Utils::DensityMatrix& densityMatrix_type2);

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
  const Eigen::MatrixXd& getCoulombMatrixType1() const;
  /**
   * @brief Getter for the Coulomb matrix
   */
  const Eigen::MatrixXd& getCoulombMatrixType2() const;

 private:
  Eigen::MatrixXd coulomb_type1_;
  Eigen::MatrixXd coulomb_type2_;
  const Utils::DensityMatrix& densityMatrix_type1_;
  const Utils::DensityMatrix& densityMatrix_type2_;

  unsigned long dim1_;
  unsigned long dim2_;
};

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_TWOTYPECOULOMBCONSTRUCTOR_H
