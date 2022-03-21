#pragma once

#include <cstdint>

namespace invertible_cdf {

/**
 * InvertibleCDF implements a reversable CDF function,
 * i.e., it is not only able to map from keys to positions
 * but also from positions to min-keys.
 */
template <class Key>
class InvertibleCDF {
 public:
  /**
   * Construct without fitting the CDF right away. Use `fit()`
   * to build the model.
   */
  InvertibleCDF() noexcept = default;

  /**
   * Construct and immediately fit CDF to data in range [begin, end)
   *
   * @tparam It iterator pointing to the data
   * @param begin iterator pointing to the first element of the sequence
   * @param end past-the-end iterator for the sequence
   */
  template <class It>
  InvertibleCDF(const It &begin, const It &end) {
    fit(begin, end);
  }

  /**
   * Fit CDF to data in range [begin, end)
   *
   * @tparam It iterator pointing to the data
   * @param begin iterator pointing to the first element of the sequence
   * @param end past-the-end iterator for the sequence
   */
  template <class It>
  void fit(const It &begin, const It &end);  // TODO(dominik): unimplemented

  /**
   * Retrieves the approximate position for a given key, i.e., computes f(x)
   *
   * @param key
   *
   * @returns position associated with `Key`
   */
  size_t pos_for_key(const Key &key);  // TODO(dominik): unimplemented

  /**
   * Retrieves the approximate key for a given position, i.e., computes
   * f^{-1}(y)
   *
   * @param pos
   *
   * @returns minimum `Key` associated with `pos`
   */
  Key key_for_pos(const size_t &pos);  // TODO(dominik): unimplemented
};
}  // namespace invertible_cdf
