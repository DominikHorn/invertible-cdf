#pragma once

#include <cstdint>
#include <vector>

#include "include/rs/builder.h"
#include "include/rs/radix_spline.h"

namespace invertible_cdf {

/**
 * InvertibleCDF implements a reversable CDF function,
 * i.e., it is not only able to map from keys to positions
 * but also from positions to min-keys.
 */
template <class Key>
class InvertibleCDF {
  /// Cdf model, i.e., f(x)
  _internal::RadixSpline<Key> rs_;

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
  void fit(const It &begin, const It &end) {
    // since we want to support arbitrarily ordered data we copy & sort
    std::vector<Key> keys(begin, end);
    std::sort(keys.begin(), keys.end());

    // use RadixSpline's builder interface for construction
    const auto min = keys.front();
    const auto max = keys.back();
    _internal::RadixSplineBuilder<Key> rsb(min, max);
    for (const auto &key : keys) rsb.AddKey(key);
    rs_ = rsb.Finalize();
  }

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

  /**
   * Equality compares two InvertibleCDF instances `a` and `b`.
   *
   * @param a
   * @param b
   *
   * @returns true iff a and b match, defined as having exactly the same
   *          internal learned cdf parameters
   */
  friend bool operator==(const InvertibleCDF<Key> &a,
                         const InvertibleCDF<Key> &b) {
    return a.rs_ == b.rs_;
  }
};
}  // namespace invertible_cdf
