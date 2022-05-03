#pragma once

#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

#include "include/builder.hpp"
#include "include/ts.hpp"

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
  plex::PLEX<Key, size_t> key_to_pos_;
  plex::PLEX<size_t, Key> pos_to_key_;

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
    plex::Builder<Key, size_t> key_to_pos_builder(
        {.x = keys.front(), .y = 0}, {.x = keys.back(), .y = keys.size() - 1},
        max_error);
    plex::Builder<size_t, Key> pos_to_key_builder(
        {.x = 0, .y = keys.front()}, {.x = keys.size() - 1, .y = keys.back()},
        max_error);

    for (size_t pos = 0; pos < keys.size(); pos++) {
      key_to_pos_builder.Add({.x = keys[pos], .y = pos});
      pos_to_key_builder.Add({.x = pos, .y = keys[pos]});
    }

    key_to_pos_ = key_to_pos_builder.Finalize();
    pos_to_key_ = pos_to_key_builder.Finalize();
  }

  /**
   * Retrieves the approximate position for a given key, i.e., computes f(x)
   *
   * @param key
   *
   * @returns bounds for positiosn associated with `key`
   */
  Bounds<size_t> pos_for_key(const Key &key) const {
    const auto search_bounds = key_to_pos_.GetSearchBound(key);
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
    const auto search_bounds = pos_to_key_.GetSearchBound(pos);
    return {.min = search_bounds.begin, .max = search_bounds.end};
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
    return a.key_to_pos_ == b.key_to_pos_ && a.pos_to_key_ == b.pos_to_key_;
  }
};
}  // namespace invertible_cdf
