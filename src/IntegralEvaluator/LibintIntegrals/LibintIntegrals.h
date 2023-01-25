/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INTEGRALEVALUATOR_LIBINTINTEGRALEVALUATOR_H
#define INTEGRALEVALUATOR_LIBINTINTEGRALEVALUATOR_H

#include <Utils/DataStructures/BasisSet.h>
#include <Utils/DataStructures/IntegralSpecifier.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Geometry/AtomCollection.h>
#include <boost/functional/hash.hpp>

namespace Scine {

namespace Utils {
class Settings;
class DensityMatrix;
} // namespace Utils
namespace Integrals {

using IntegralEvaluatorMap =
    std::unordered_map<Utils::Integrals::ReturnKey, Eigen::MatrixXd, boost::hash<Utils::Integrals::ReturnKey>>;

class LibintIntegrals {
 public:
  static constexpr const char* model = "libint_integrals";

  LibintIntegrals();

  ~LibintIntegrals();

  /**
   * @brief Evaluate integrals between two basis sets.
   * Operator and other infos are contained in the specifier object.
   * @return unordered_map containing matrices of type Eigen::MatrixXd
   */
  static auto evaluate(const Utils::Integrals::IntegralSpecifier& specifier, const Utils::Integrals::BasisSet& basis1,
                       const Utils::Integrals::BasisSet& basis2) -> IntegralEvaluatorMap;
  /**
   * @brief Simple initialization routine of a basis set object, that uses libint functionality.
   * Must be passed the name of the basis set and the atoms object.
   * The libint routine searches in the basis set file with the given `name` for the basis set specifics of the
   * atoms that are contained in the `atoms` object.
   * @param name
   * @param atoms
   * @return Scine BasisSet object
   */
  auto initializeBasisSet(const std::string& name, const Utils::AtomCollection& atoms, bool generateShellPairs = true)
      -> Utils::Integrals::BasisSet;
  /**
   * @brief Initializes a BasisSet located at the atoms' positions.
   *        Used for initializing different basis set for all atoms.
   * @param name of the basis set, e.g., def2-TZVP, def2-SVP, ...
   * @param atoms - contains positions and charges of the atoms
   * @return BasisSet object
   */
  auto initializeBasisSet(const std::unordered_map<Utils::ElementType, std::string>& names,
                          const Utils::AtomCollection& atoms, bool generateShellPairs = true) -> Utils::Integrals::BasisSet;

  /**
   * @brief Computes non-negligible shell pair list; shells \c i and \c j form a
   * non-negligible pair if they share a center or the Frobenius norm of their overlap is
   * greater than threshold.
   * @param performOverlapPrescreening This flags whether a screening of pairs has to be done based on the
   *                                   pseudo density \f$ (\mu \nu |\f$.
   * @param threshold Overlap threshold for computing relevant shell-pairs. Defaults to 1e-12.
   * @param calculateCauchySchwarzFactor Decides whether the maximal Cauchy-Schwarz factor for the shell-pair
   *                                     has to be calculated. Needed to do a Cauchy-Schwartz ERI prescreening.
   */
  static auto generateShellPairs(Utils::Integrals::BasisSet& basis, bool performOverlapPrescreening = true,
                                 double threshold = 1e-12, bool calculateCauchySchwarzFactor = true) -> void;

  /**
   * @brief Evaluates the pre-BO contribution to the Fock matrix for different particle types.
   * Operator and other infos are contained in the specifier object.
   * @param specifier
   * @param basis1
   * @param basis2
   * @param dm1
   * @param dm2
   * @return (Matrix corresponding to dm1 contraction, Matrix corresponding to dm2 contraction)
   */
  static auto evaluateTwoBodyDirectPreBo(const Utils::Integrals::IntegralSpecifier& specifier,
                                         const Utils::Integrals::BasisSet& basis1,
                                         const Utils::Integrals::BasisSet& basis2, const Utils::DensityMatrix& dm1,
                                         const Utils::DensityMatrix& dm2) -> std::pair<Eigen::MatrixXd, Eigen::MatrixXd>;

  /**
   * @brief Evaluates the BO contribution to the Fock matrix for different particle types.
   * Cauchy-Schwarz pre-screening is performed if basis1==basis2
   * Operator and other infos are contained in the specifier object.
   * @param specifier
   * @param basis1
   * @param basis2
   * @param dm1
   * @param prescreeningThreshold
   * @return J, K matrices.
   */
  static auto evaluateTwoBodyDirectBo(const Utils::Integrals::IntegralSpecifier& specifier,
                                      const Utils::Integrals::BasisSet& basis1, const Utils::Integrals::BasisSet& basis2,
                                      const Utils::DensityMatrix& dm1, double prescreeningThreshold)
      -> std::pair<Utils::SpinAdaptedMatrix, Utils::SpinAdaptedMatrix>;
  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Utils::Settings& settings();
  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Utils::Settings& settings() const;

  static std::string name();

 private:
  bool _verbose = true;
  std::unique_ptr<Utils::Settings> _settings;

  /**
   * @brief Evaluate one-body integrals between two basis sets.
   * Operator and other infos are contained in the specifier object.
   * @return unordered_map containing matrices of type Eigen::MatrixXd
   */
  static auto evaluateOneBody(const Utils::Integrals::IntegralSpecifier& specifier, const Utils::Integrals::BasisSet& basis1,
                              const Utils::Integrals::BasisSet& basis2) -> IntegralEvaluatorMap;
  /**
   * @brief Evaluate two-body integrals between two basis sets.
   * Operator and other infos are contained in the specifier object.
   * @return unordered_map containing matrices of type Eigen::MatrixXd
   */
  static auto evaluateTwoBody(const Utils::Integrals::IntegralSpecifier& specifier, const Utils::Integrals::BasisSet& basis1,
                              const Utils::Integrals::BasisSet& basis2) -> IntegralEvaluatorMap;
};

} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_LIBINTINTEGRALEVALUATOR_H
