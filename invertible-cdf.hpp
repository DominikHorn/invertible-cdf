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
template <class Key, size_t max_error = 16>
class InvertibleCDF {
  /// Cdf model, i.e., f(x)
  _internal::RadixSpline<Key> rs_;

 public:
  /**
   * Construct without fitting the CDF right away. Use `train()`
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
    train(begin, end);
  }

  /**
   * Fit CDF to data in range [begin, end)
   *
   * @tparam It iterator pointing to the data
   * @param begin iterator pointing to the first element of the sequence
   * @param end past-the-end iterator for the sequence
   */
  template <class It>
  void train(const It &begin, const It &end) {
    // since we want to support arbitrarily ordered data we copy & sort
    std::vector<Key> keys(begin, end);
    std::sort(keys.begin(), keys.end());

    // use RadixSpline's builder interface for construction
    const auto min = keys.front();
    const auto max = keys.back();
    _internal::RadixSplineBuilder<Key> rsb(min, max, 1, max_error);
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
   * @returns bounds for keys associated with `pos`, if `pos` in range indexed
   * by cdf. Otherwise `{ .min = Key::max(), .max = Key::max() }`
   */
  Bounds<Key> key_for_pos(const size_t &pos) const {
    using KeyLims = std::numeric_limits<Key>;
    using SegmentIter = typename decltype(rs_.spline_points_)::const_iterator;

    // rename to enhance below code's readability
    const auto &spline = rs_.spline_points_;

    // min and max pos with error bound taken into account. Our true
    // key could map to any position in this range -> we need to
    // extrapolate from this raneg
    const auto min_pos = pos - std::min(max_error, pos);
    const auto max_pos =
        pos + std::min(rs_.pos_table_.size() - 1 - pos, max_error);

    // find first segment with seg.y <= pos and return it's key (x)
    const auto min_seg = spline.begin() + rs_.pos_table_[min_pos];
    const auto max_seg = spline.begin() + rs_.pos_table_[max_pos];

    assert(min_seg >= spline.begin());
    assert(min_seg < spline.end());
    assert(max_seg >= spline.begin());
    assert(max_seg < spline.end());
    assert(min_seg <= max_seg);

    // found segment for our pos
    const auto inverted_fma = [](const SegmentIter &up, const SegmentIter &down,
                                 const size_t p) {
      const double slope = (up->y - down->y) / (up->x - down->x);
      const double intercept = down->y;

      // y = ax + b <-> x = (y-b)/a
      // NOTE: x on slope is in local space and not global space, hence add
      // prev_segment->x
      return down->x + (static_cast<double>(p) - intercept) / slope;
    };

    const double pred_min =
        std::floor(inverted_fma(min_seg + 1, min_seg, min_pos));

    // key was never inserted during train(), i.e., is 'beyond' our spline
    if (max_seg + 1 >= spline.end()) {
      return {.min = static_cast<Key>(pred_min), .max = KeyLims::max()};
    }

    const double pred_max =
        std::ceil(inverted_fma(max_seg + 1, max_seg, max_pos));

    return {.min = static_cast<Key>(pred_min),
            .max = static_cast<Key>(pred_max)};
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
  friend bool operator==(const InvertibleCDF<Key, max_error> &a,
                         const InvertibleCDF<Key, max_error> &b) {
    return a.rs_ == b.rs_;
  }
};
}  // namespace invertible_cdf
