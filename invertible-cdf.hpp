#pragma once

#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

#include "include/rs/builder.h"
#include "include/rs/radix_spline.h"

namespace invertible_cdf {

/**
 * Bounds represent a [min, max) interval in a certain domain
 *
 * @tparam T domain, e.g., std::size_t
 */
template <class T>
struct Bounds {
  T min, max;
};

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
   * @returns bounds for positiosn associated with `key`
   */
  Bounds<size_t> pos_for_key(const Key &key) const {
    const auto search_bounds = rs_.GetSearchBound(key);
    return {.min = search_bounds.begin, .max = search_bounds.end};
  }

  /**
   * Retrieves the approximate key for a given position, i.e., computes
   * f^{-1}(y)
   *
   * @param pos
   *
   * @returns bounds for keys associated with `pos`, if `\pos` in range indexed
   * by cdf. Otherwise `{ .min = Key::max(), .max = Key::max() }`
   */
  Bounds<Key> key_for_pos(const size_t &pos) const {
    using KeyLims = std::numeric_limits<Key>;

    // rename to enhance below code's readability
    const auto &spline = rs_.spline_points_;

    // find first segment with seg.y <= pos and return it's key (x)
    for (auto curr_segment = spline.begin() + 1, prev_segment = spline.begin();
         curr_segment < spline.end(); prev_segment++, curr_segment++) {
      // return as soon as we find the position's lower bound
      if (curr_segment->y >= pos) {
        return {.min = prev_segment->x, .max = curr_segment->x};
      }
    }

    return {.min = spline.back().x, .max = KeyLims::max()};
  }

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
