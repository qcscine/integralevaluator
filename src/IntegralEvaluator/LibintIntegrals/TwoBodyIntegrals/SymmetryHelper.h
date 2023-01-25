/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INTEGRALEVALUATOR_SYMMETRYHELPER_H
#define INTEGRALEVALUATOR_SYMMETRYHELPER_H

#include <array>
#include <vector>

#define UNUSED(expr) \
  do {               \
    (void)(expr);    \
  } while (0)

namespace Scine {
namespace Integrals {
namespace TwoBody {

template<int rank>
using IndexType = std::array<int, rank>;

enum class IntegralSymmetry { onefold, twofold, fourfold, eightfold };

/**
 * @brief Function to get the degeneracy of a index quartet with some symmetry.
 * @tparam symmetry The symmetry. Can be onefold (no symmetry), twofold, fourfold or eightfold.
 * @param basisFunction1
 * @param basisFunction2
 * @param basisFunction3
 * @param basisFunction4
 * @return A double giving the degeneracy of the index quartet in the given symmetry.
 */
template<IntegralSymmetry symmetry>
double getDegeneracy(int basisFunction1, int basisFunction2, int basisFunction3, int basisFunction4);

/**
 * @brief Maps an unsymmetrized two-electron index to the right element in a symmetrized ERI map.
 * @tparam symmetry The symmetry of the map.
 * @param index The unsymmetrized index quartet.
 * @return The index quartet pointing to the corresponding index.
 */
template<IntegralSymmetry symmetry>
IndexType<4> getMappedIndex(IndexType<4> index);

/**
 * @brief Returns a pack of indices related by symmetry to the input one.
 * @tparam symmetry The Eri symmetry.
 * @param index The symmetrized index.
 * @return A vector containing the indices related by symmetry to the input one.
 */
template<IntegralSymmetry symmetry>
inline std::vector<IndexType<4>> getSymmetricIndices(IndexType<4> index);

/**
 * @brief one-to-one mapping of the symmetrized index to itself.
 * @param index the symmetrized index.
 * @return A vector containing the same index.
 * Since there is no symmetry in onefold symmetrized ERIs,
 * the return value is the same as the input one.
 * Example:
 * Symmetrized index: 3210
 * Indices:           3210
 */
template<>
inline std::vector<IndexType<4>> getSymmetricIndices<IntegralSymmetry::onefold>(IndexType<4> index) {
  return {index};
}

/**
 * @brief one-to-many mapping of the symmetrized index to the actual ones.
 * @param index the symmetrized index.
 * @return A vector containing all the indices that are symmetric to the input one.
 * Example:
 * Symmetrized index: 3210
 * Indices:           3210 1032
 * Symmetrized index: 2100
 * Indices:           2100 0021
 */
template<>
inline std::vector<IndexType<4>> getSymmetricIndices<IntegralSymmetry::twofold>(IndexType<4> index) {
  if (index[0] != index[2] || index[1] != index[3]) {
    return {index, {index[2], index[3], index[0], index[1]}};
  }
  return {index};
}

/**
 * @brief one-to-many mapping of the symmetrized index to the actual ones.
 * @param index the symmetrized index.
 * @return A vector containing all the indices that are symmetric to the input one.
 * Example:
 * Symmetrized  index: 3210
 * Indices:            3210 3201 2310 2301
 * Symmetrized index:  2100
 * Indices:            2100 1200
 */
template<>
inline std::vector<IndexType<4>> getSymmetricIndices<IntegralSymmetry::fourfold>(IndexType<4> index) {
  std::vector<IndexType<4>> result;
  result.reserve(4);
  result.push_back(index);
  if (index[0] != index[1]) {
    result.push_back({index[1], index[0], index[2], index[3]});
  }
  if (index[2] != index[3]) {
    result.push_back({index[0], index[1], index[3], index[2]});
  }
  if (index[0] != index[1] && index[2] != index[3]) {
    result.push_back({index[1], index[0], index[3], index[2]});
  }
  return result;
}
/**
 * @brief one-to-many mapping of the symmetrized index to the actual ones.
 * @param index the symmetrized index.
 * @return A vector containing all the indices that are symmetric to the input one.
 * Example:
 * Symmetrized index: 3210
 * Indices:           3210 3201 2310 2301 1032 0132 1023 0123
 * Symmetrized index: 3232
 * Indices:           3232 3223 2332 2323
 * Symmetrized index: 2100
 * Indices:           2100 1200 0021 0012
 * Symmetrized index: 1100
 * Indices:           1100 0011
 * Symmetrized index: 0000
 * Indices:           0000
 */
template<>
inline std::vector<IndexType<4>> getSymmetricIndices<IntegralSymmetry::eightfold>(IndexType<4> index) {
  std::vector<IndexType<4>> result;
  result.reserve(8);
  result.push_back(index);

  if (!(index[0] == index[1] && index[0] == index[2] && index[0] == index[3])) {
    if (index[0] == index[1] && index[2] == index[3]) { // 1100
      result.emplace_back(IndexType<4>{index[2], index[2], index[0], index[0]});
    }
    else if (index[0] == index[1] && index[2] != index[3]) { // 1110
      result.emplace_back(IndexType<4>{index[0], index[0], index[3], index[2]});
      result.emplace_back(IndexType<4>{index[2], index[3], index[0], index[0]});
      result.emplace_back(IndexType<4>{index[3], index[2], index[0], index[0]});
    }
    else if (index[2] == index[3] && index[0] != index[1]) { // 2100
      result.emplace_back(IndexType<4>{index[1], index[0], index[2], index[2]});
      result.emplace_back(IndexType<4>{index[2], index[2], index[0], index[1]});
      result.emplace_back(IndexType<4>{index[2], index[2], index[1], index[0]});
    }
    else if (index[0] == index[2] && index[1] == index[3] && index[0] != index[1]) { // 3232
      result.emplace_back(IndexType<4>{index[1], index[0], index[2], index[3]});
      result.emplace_back(IndexType<4>{index[0], index[1], index[3], index[2]});
      result.emplace_back(IndexType<4>{index[1], index[0], index[3], index[2]});
    }
    else if (index[0] != index[1] && index[2] != index[3] && index[0] != index[3]) { // 3210
      result.emplace_back(IndexType<4>{index[1], index[0], index[2], index[3]});
      result.emplace_back(IndexType<4>{index[0], index[1], index[3], index[2]});
      result.emplace_back(IndexType<4>{index[1], index[0], index[3], index[2]});
      result.emplace_back(IndexType<4>{index[2], index[3], index[0], index[1]});
      result.emplace_back(IndexType<4>{index[2], index[3], index[1], index[0]});
      result.emplace_back(IndexType<4>{index[3], index[2], index[0], index[1]});
      result.emplace_back(IndexType<4>{index[3], index[2], index[1], index[0]});
    }
  }
  return result;
}
template<>
inline double getDegeneracy<IntegralSymmetry::onefold>(int basisFunction1, int basisFunction2, int basisFunction3,
                                                       int basisFunction4) {
  UNUSED(basisFunction1);
  UNUSED(basisFunction2);
  UNUSED(basisFunction3);
  UNUSED(basisFunction4);
  return 1.0;
}
template<>
inline double getDegeneracy<IntegralSymmetry::twofold>(int basisFunction1, int basisFunction2, int basisFunction3,
                                                       int basisFunction4) {
  return (basisFunction1 == basisFunction3) ? (basisFunction2 == basisFunction4 ? 1.0 : 2.0) : 2.0;
}
template<>
inline double getDegeneracy<IntegralSymmetry::fourfold>(int basisFunction1, int basisFunction2, int basisFunction3,
                                                        int basisFunction4) {
  auto degeneracy_12 = (basisFunction1 == basisFunction2) ? 1.0 : 2.0;
  auto degeneracy_34 = (basisFunction3 == basisFunction4) ? 1.0 : 2.0;
  return degeneracy_12 * degeneracy_34;
}
template<>
inline double getDegeneracy<IntegralSymmetry::eightfold>(int basisFunction1, int basisFunction2, int basisFunction3,
                                                         int basisFunction4) {
  auto degeneracy_fourfold =
      getDegeneracy<IntegralSymmetry::fourfold>(basisFunction1, basisFunction2, basisFunction3, basisFunction4);
  auto degeneracy_twofold =
      getDegeneracy<IntegralSymmetry::twofold>(basisFunction1, basisFunction2, basisFunction3, basisFunction4);
  return degeneracy_fourfold * degeneracy_twofold;
}

template<>
inline IndexType<4> getMappedIndex<IntegralSymmetry::onefold>(IndexType<4> index) {
  return index;
}

template<>
inline IndexType<4> getMappedIndex<IntegralSymmetry::twofold>(IndexType<4> index) {
  if (index[2] > index[0] || (index[2] == index[0] && index[3] > index[1])) {
    std::swap(index[2], index[0]);
    std::swap(index[3], index[1]);
  }
  return index;
}

template<>
inline IndexType<4> getMappedIndex<IntegralSymmetry::fourfold>(IndexType<4> index) {
  if (index[3] > index[2])
    std::swap(index[3], index[2]);
  if (index[1] > index[0])
    std::swap(index[1], index[0]);
  return index;
}

template<>
inline IndexType<4> getMappedIndex<IntegralSymmetry::eightfold>(IndexType<4> index) {
  if (index[3] > index[2])
    std::swap(index[3], index[2]);
  if (index[1] > index[0])
    std::swap(index[1], index[0]);
  if (index[2] > index[0] || (index[2] == index[0] && index[3] > index[1])) {
    std::swap(index[2], index[0]);
    std::swap(index[3], index[1]);
  }
  return index;
}

} // namespace TwoBody
} // namespace Integrals
} // namespace Scine

#endif // INTEGRALEVALUATOR_SYMMETRYHELPER_H
