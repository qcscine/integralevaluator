/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_LIBINT_H
#define INTEGRALEVALUATOR_LIBINT_H

/* Include Std and External Headers */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <libint2.hpp>
#pragma GCC diagnostic pop

#include <Utils/DataStructures/BasisSet.h>

namespace Scine {
namespace Integrals {
/**
 * @class Libint @file Libint.h
 * @brief Singleton to initialize and wrap up Libint.
 * This class assumes that the maximal number of primitives per contracted Gaussian
 * shell is 20.
 */
class Libint {
 public:
  /**
   * Deleted functions neede for the singleton pattern.
   * @{
   */
  Libint(const Libint&) = delete;
  Libint& operator=(const Libint&) = delete;
  //! @}
 private:
  /**
   * @brief Private constructor, part of the singleton pattern.
   */
  Libint();

 public:
  /**
   * @brief Getter for an instance of the Libint singleton.
   * The private constructor must initialize the Libint environment.
   */
  static Libint& getInstance() {
    static Libint libint;
    return libint;
  }

  /**
   * @brief Function returning the max number of threads in use.
   * @return the value of omp_get_max_threads() at initialization time of this class.
   */
  static int getMaxNumberThreads() {
    return nThreads_;
  }
  /**
   * @brief Destructor. Must finalize Libint.
   */
  ~Libint();

  static auto getEngine(const Scine::Utils::Integrals::BasisSet& basis, const libint2::Operator& op,
                        const int& derivOrder = 0) -> libint2::Engine;

  static auto getEngine(const Scine::Utils::Integrals::BasisSet& basis1, const Scine::Utils::Integrals::BasisSet& basis2,
                        const libint2::Operator& op, const int& derivOrder) -> libint2::Engine;

  static auto getEngine(const libint2::Operator& op, const libint2::Shell& shell1, const libint2::Shell& shell2,
                        const libint2::Shell& shell3, const libint2::Shell& shell4) -> libint2::Engine;

 private:
  // Dependent on the basis.
  static constexpr int maxPrimitiveNumber = 20;
  // Dependent on how the libint2 library was generated. Use Libint2 macro
  static constexpr int maxAngularMomentum = LIBINT2_MAX_AM;
  // Saves the number of OMP threads to use
  static int nThreads_;
};

} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_LIBINT_H
