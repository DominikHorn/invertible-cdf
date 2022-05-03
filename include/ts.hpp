#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "cht/cht.hpp"
#include "common.hpp"

namespace invertible_cdf::plex {
// Approximates a cumulative distribution function (CDF) using spline
// interpolation.
template <class X, class Y>
class PLEX {
  using Coord = Coord<X, Y>;

 public:
  PLEX() = default;

  PLEX(const Coord& min, const Coord& max, size_t num_keys,
       size_t spline_max_error, cht::CompactHistTree<X> cht,
       std::vector<Coord> spline_points)
      : min_(min),
        max_(max),
        num_coords_(num_keys),
        spline_max_error_(spline_max_error),
        spline_points_(std::move(spline_points)),
        cht_(std::move(cht)) {}

  double GetEstimatedY(const X x) const {
    // Truncate to data boundaries.
    if (x <= min_.x) return min_.y;
    if (x >= max_.x) return max_.y;

    // Find spline segment with `key` ∈ (spline[index - 1], spline[index]].
    const size_t index = GetSplineSegment(x);
    const Coord down = spline_points_[index - 1];
    const Coord up = spline_points_[index];

    // Compute slope.
    assert(up.x >= down.x);
    assert(up.y >= down.y);
    const double x_diff = up.x - down.x;
    const double y_diff = up.y - down.y;
    const double slope = y_diff / x_diff;

    // Interpolate.
    const double key_diff = static_cast<double>(x) - static_cast<double>(down.x);
    return std::fma(key_diff, slope, down.y);
  }

  // Returns a search bound [begin, end) around the estimated position.
  SearchBound GetSearchBound(const X x) const {
    const size_t estimate = GetEstimatedY(x);
    const size_t begin = estimate - std::min(static_cast<std::uint64_t>(estimate), static_cast<std::uint64_t>(spline_max_error_));
    // `end` is exclusive.
    const size_t end = estimate + std::min(static_cast<std::uint64_t>(max_.y - estimate), static_cast<std::uint64_t>(spline_max_error_));
    return SearchBound{begin, end};
  }

  // Returns the size in bytes.
  size_t GetSize() const {
    return sizeof(*this) + cht_.GetSize() +
           spline_points_.size() * sizeof(Coord);
  }

  friend bool operator==(const PLEX<X, Y>& a, const PLEX<X, Y>& b) {
    return a.spline_points_ == b.spline_points_;
  }

 private:
  // Returns the index of the spline point that marks the end of the spline
  // segment that contains the `key`: `key` ∈ (spline[index - 1], spline[index]]
  size_t GetSplineSegment(const X key) const {
    // Narrow search range using CHT.
    const auto range = cht_.GetSearchBound(key);

    // Linear search?
    if (range.end - range.begin < 32) {
      // Do linear search over narrowed range.
      uint32_t current = range.begin;
      while (spline_points_[current].x < key) ++current;
      return current;
    }

    // Do binary search over narrowed range.
    const auto lb = std::lower_bound(
        spline_points_.begin() + range.begin,
        spline_points_.begin() + range.end, key,
        [](const Coord& coord, const X key) { return coord.x < key; });
    return std::distance(spline_points_.begin(), lb);
  }

  Coord min_;
  Coord max_;
  size_t num_coords_;
  size_t spline_max_error_;

  std::vector<Coord> spline_points_;
  cht::CompactHistTree<X> cht_;
};

}  // namespace invertible_cdf::plex
