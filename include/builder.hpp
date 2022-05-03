#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <optional>

#include "cht/builder.hpp"
#include "cht/cht.hpp"
#include "common.hpp"
#include "ts.hpp"

namespace invertible_cdf::plex {
// Builds a `TrieSpline`.
template <class X, class Y>
class Builder {
  using Coord = Coord<X, Y>;

 public:
  Builder(const Coord& min, const Coord& max, Y spline_max_error)
      : min_(min),
        max_(max),
        spline_max_error_(spline_max_error),
        curr_num_coords_(0),
        curr_num_distinct_coords_(0),
        prev_coord_(min),
        chtb_(min.x, max.x) {}

  void Add(const Coord& coord) {
    assert(coord.x >= min_.x && coord.x <= max_.x);
    assert(coord.y >= min_.y && coord.y <= max_.y);

    // Inputs need to be monotonically increasing.
    assert(coord.x >= prev_coord_.x);
    assert(coord.y >= prev_coord_.y);

    PossiblyAddCoordToSpline(coord);

    ++curr_num_coords_;
    prev_coord_ = coord;
  }

  // Finalizes the construction and returns a read-only `TrieSpline`.
  PLEX<X, Y> Finalize() {
    // Last coord needs to be equal to `max_`.
    assert(curr_num_coords_ == 0 || prev_coord_ == max_);

    // Ensure that `prev_key_` (== `max_key_`) is last key on spline.
    if (curr_num_coords_ > 0 && spline_points_.back() != prev_coord_) {
      AddCoordToSpline(prev_coord_);
    }

    // Find tuning.
    std::vector<Statistics> statistics;
    ComputeStatistics(statistics);
    auto tuning = InferTuning(statistics);

    // Finalize CHT
    auto cht_ = chtb_.Finalize(tuning.numBins, tuning.treeMaxError);

    // And return the read-only instance
    return PLEX<X, Y>(min_, max_, curr_num_coords_, spline_max_error_,
                      std::move(cht_), std::move(spline_points_));
  }

 private:
  using Interval = std::pair<unsigned, unsigned>;
  using Statistics = Statistics;

  template <class T>
  static unsigned ComputeLog(T n, bool round = false) {
    assert(n);
    auto lz = 0;
    if constexpr (sizeof(T) < 8) {
      lz = __builtin_clz(n);
    } else {
      lz = __builtin_clzl(n);
    }

    return (8 * sizeof(T) - 1) - lz + (round ? ((n & (n - 1)) != 0) : 0);
  }

  template <class T>
  static unsigned ComputeLcp(T x, T y) {
    if constexpr (sizeof(T) * 8 < 64) return __builtin_clz(x ^ y);
    return __builtin_clzl(x ^ y);
  }

  template <class T>
  static size_t GetNumShiftBits(T diff, size_t num_radix_bits) {
    const uint32_t clzl = __builtin_clzl(diff);
    if ((8 * sizeof(T) - clzl) < num_radix_bits) return 0;
    return 8 * sizeof(T) - num_radix_bits - clzl;
  }

  // `int(ceil(log_2(distance)))`. I
  static size_t ComputeCost(uint32_t value) {
    assert(value);
    if (value == 1) return 1;
    return 31 - __builtin_clz(value) + ((value & (value - 1)) != 0);
  }

  void AddCoordToSpline(const Coord& coord) {
    spline_points_.push_back({coord.x, coord.y});
    AddCoordToCHT(coord);
  }

  enum Orientation { Collinear, CW, CCW };
  static constexpr double precision = std::numeric_limits<double>::epsilon();

  static Orientation ComputeOrientation(const double dx1, const double dy1,
                                        const double dx2, const double dy2) {
    const double expr = std::fma(dy1, dx2, -std::fma(dy2, dx1, 0));
    if (expr > precision) {
      return Orientation::CW;
    } else if (expr < -precision) {
      return Orientation::CCW;
    }
    return Orientation::Collinear;
  };

  // Implementation is based on `GreedySplineCorridor` from:
  // T. Neumann and S. Michel. Smooth interpolating histograms with error
  // guarantees. [BNCOD'08]
  void PossiblyAddCoordToSpline(const Coord& coord) {
    if (curr_num_coords_ == 0) {
      // Add first CDF point to spline.
      AddCoordToSpline(coord);
      ++curr_num_distinct_coords_;
      prev_coord_ = coord;
      return;
    }

    if (coord.x == prev_coord_.x) {
      // No new CDF point if x didn't change.
      return;
    }

    // New CDF point.
    ++curr_num_distinct_coords_;

    // Compute current `upper_y` and `lower_y`.
    const auto lower_y = coord.y - std::min(coord.y, spline_max_error_);
    const auto upper_y = coord.y + spline_max_error_;

    if (curr_num_distinct_coords_ == 2) {
      // Initialize `upper_limit_` and `lower_limit_` using the second CDF
      // point.
      lower_limit_ = {.x = coord.x, .y = lower_y};
      upper_limit_ = {.x = coord.x, .y = upper_y};
      prev_coord_ = coord;
      return;
    }

    // `B` in algorithm.
    const Coord& last = spline_points_.back();

    assert(upper_limit_.x >= last.x);
    assert(upper_limit_.y >= last.y);
    assert(lower_limit_.x >= last.x);
    assert(coord.x >= last.x);
    assert(coord.y >= last.y);

    // Compute differences.
    const double upper_limit_x_diff =
        static_cast<double>(upper_limit_.x) - static_cast<double>(last.x);
    const double lower_limit_x_diff =
        static_cast<double>(lower_limit_.x) - static_cast<double>(last.x);
    const double x_diff = coord.x - last.x;

    const double upper_limit_y_diff =
        static_cast<double>(upper_limit_.y) - static_cast<double>(last.y);
    const double lower_limit_y_diff =
        static_cast<double>(lower_limit_.y) - static_cast<double>(last.y);
    const double y_diff = coord.y - last.y;

    // `prev_coord_` is the previous point on the CDF and the next candidate to
    // be added to the spline. Hence, it should be different from the `last`
    // point on the spline.
    assert(prev_coord_ != last);

    // Do we cut the error corridor?
    if ((ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff,
                            y_diff) != Orientation::CW) ||
        (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff,
                            y_diff) != Orientation::CCW)) {
      // Add previous CDF point to spline.
      AddCoordToSpline(prev_coord_);

      // Update limits.
      lower_limit_ = {.x = coord.x, .y = lower_y};
      upper_limit_ = {.x = coord.x, .y = upper_y};
    } else {
      assert(upper_y >= last.y);

      const double upper_y_diff = upper_y - last.y;
      if (ComputeOrientation(upper_limit_x_diff, upper_limit_y_diff, x_diff,
                             upper_y_diff) == Orientation::CW) {
        upper_limit_ = {.x = coord.x, .y = upper_y};
      }

      const double lower_y_diff = lower_y - last.y;
      if (ComputeOrientation(lower_limit_x_diff, lower_limit_y_diff, x_diff,
                             lower_y_diff) == Orientation::CCW) {
        lower_limit_ = {.x = coord.x, .y = lower_y};
      }
    }

    prev_coord_ = coord;
  }

  void AddCoordToCHT(const Coord& coord) { chtb_.AddKey(coord.x); }

  void ComputeRadixTableStatistics(std::vector<Statistics>& statistics) {
    static constexpr unsigned maxNumRadixBits = 30;

    // Init the radix analyzer.
    std::vector<RadixConfig> radixAnalyzer;
    radixAnalyzer.resize(1 + maxNumRadixBits);
    for (unsigned index = 1; index <= maxNumRadixBits; ++index) {
      radixAnalyzer[index].shiftBits = GetNumShiftBits(max_.x - min_.x, index);
      radixAnalyzer[index].prevPrefix = 0;
      radixAnalyzer[index].prevSplineIndex = 0;
      radixAnalyzer[index].cost = 0;
    }

    // And compute the costs.
    for (unsigned splineIndex = 1, limit = spline_points_.size();
         splineIndex != limit; ++splineIndex) {
      for (unsigned radix = 1; radix <= maxNumRadixBits; ++radix) {
        const X currPrefix = (spline_points_[splineIndex].x - min_.x) >>
                             radixAnalyzer[radix].shiftBits;

        // New prefix?
        if (currPrefix != radixAnalyzer[radix].prevPrefix) {
          // Then compute statistics.
          assert(splineIndex);
          const unsigned prevSplineIndex = radixAnalyzer[radix].prevSplineIndex;
          const size_t numDataKeys =
              spline_points_[splineIndex].y - spline_points_[prevSplineIndex].y;
          const size_t numSplineKeys = splineIndex - prevSplineIndex;
          assert(numSplineKeys);

          // Update cost.
          radixAnalyzer[radix].cost += numDataKeys * ComputeCost(numSplineKeys);

          // Update the parameters.
          radixAnalyzer[radix].prevPrefix = currPrefix;
          radixAnalyzer[radix].prevSplineIndex = splineIndex;
        }
      }
    }

    // Finalize the costs.
    for (unsigned radix = 1; radix <= maxNumRadixBits; ++radix) {
      // Compute statistics.
      const unsigned prevSplineIndex = radixAnalyzer[radix].prevSplineIndex;
      const size_t numDataKeys =
          spline_points_.back().y - spline_points_[prevSplineIndex].y;
      const size_t numSplineKeys = spline_points_.size() - prevSplineIndex;
      assert(numSplineKeys);

      // Update the cost.
      radixAnalyzer[radix].cost += numDataKeys * ComputeCost(numSplineKeys);

      // Normalize the cost.
      radixAnalyzer[radix].cost /= spline_points_.back().y;

      // And save them into `statistics`.
      statistics.emplace_back(
          1u << radix, std::numeric_limits<unsigned>::max(),
          radixAnalyzer[radix].cost,
          static_cast<size_t>(
              ((max_.x - min_.x) >> radixAnalyzer[radix].shiftBits) + 2) *
              sizeof(unsigned));
    }
  }

  void ComputeCHTStatistics(std::vector<Statistics>& statistics) {
    // Compute the necessary amount of bits we need.
    const unsigned lg = ComputeLog(max_.x - min_.x, true);
    const unsigned alreadyCommon = (sizeof(X) << 3) - lg;

    // Compute the longest common prefix.
    const auto ExtractLCP = [&](unsigned index) -> unsigned {
      const auto lcp = ComputeLcp(spline_points_[index].x - min_.x,
                                  spline_points_[index - 1].x - min_.x);
      assert(lcp >= alreadyCommon);
      return lcp - alreadyCommon;  // __builtin_clzl((spline_points_[index].x -
                                   // min_key_) ^ (spline_points_[index - 1].x -
                                   // min_key_)) - alreadyCommon;
    };

    // Fill the lcp-array with lcp[i] := lcp(key[i], key[i - 1]).
    std::vector<unsigned> lcp(spline_points_.size());
    std::vector<unsigned> counters(1 + lg);
    lcp[0] = std::numeric_limits<unsigned>::max();
    for (unsigned index = 1, limit = spline_points_.size(); index != limit;
         ++index) {
      lcp[index] = ExtractLCP(index);
      counters[lcp[index]]++;
    }

    // Sort the lcp-array.
    std::vector<unsigned> offsets(1 + lg);
    for (unsigned bit = 1; bit <= lg; ++bit) {
      offsets[bit] = offsets[bit - 1] + counters[bit - 1];
    }
    std::vector<unsigned> sorted(spline_points_.size() - 1);
    for (unsigned index = 1, limit = spline_points_.size(); index != limit;
         ++index) {
      sorted[offsets[lcp[index]]++] = index;
    }

    // Init the histogram. We only store the current row and the previous row,
    // as we do not need all others.
    std::pair<unsigned, std::vector<Interval>> histogram[2];
    histogram[0].first = histogram[1].first = 0;
    histogram[0].second.resize(spline_points_.size());
    histogram[1].second.resize(spline_points_.size());
    unsigned side = 0;

    // Add a new interval.
    const auto AddNewInterval = [&](Interval interval) -> void {
      // Empty interval?
      if (interval.first == interval.second) return;
      histogram[side].second[histogram[side].first++] = interval;
    };

    // Init the first level of the histogram.
    AddNewInterval({1, spline_points_.size()});
    unsigned ptrInSorted = 0;

    const auto AnalyzeInterval = [&](unsigned level,
                                     Interval interval) -> void {
      // Check whether `pos` lies inside `interval`.
      const auto isInside = [&](unsigned pos) -> bool {
        return (interval.first <= pos) && (pos < interval.second);
      };

      // Already fully-consumed?
      if (ptrInSorted == sorted.size()) return;

      // Could the next `lcp` be an interval-breaker?
      auto leftSide = interval.first;
      while ((ptrInSorted != sorted.size()) &&
             (lcp[sorted[ptrInSorted]] < level)) {
        // Not inside?
        if (!isInside(sorted[ptrInSorted])) break;

        // Otherwise, split the interval and create a new one.
        AddNewInterval({leftSide, sorted[ptrInSorted]});

        // Reset the left side.
        leftSide = sorted[ptrInSorted] + 1;
        ++ptrInSorted;
      }

      // And add the remaining interval. It could be as the original one if
      // there were no interval-breakers.
      AddNewInterval({leftSide, interval.second});
    };

    // Init the statistics.
    static constexpr unsigned maxNumPossibleBins = 20;
    static constexpr unsigned maxPossibleTreeError = 1u << 10;
    const unsigned numPossibleBins = std::min(maxNumPossibleBins, lg);
    std::vector<unsigned> possibleNumBins(numPossibleBins);
    std::vector<std::vector<std::pair<uint64_t, unsigned>>> matrix(
        numPossibleBins);

    // Set the possible number of bins and init the matrix which represents the
    // prefix sums.
    for (unsigned index = 0; index != numPossibleBins; ++index) {
      possibleNumBins[index] = (1u << (index + 1));
      matrix[index].assign(1 + maxPossibleTreeError, {0, 0});
    }

    // Consume the `level`th row of the histogram.
    const auto ConsumeLevel = [&](unsigned level) -> void {
      for (unsigned index = 0; index != numPossibleBins; ++index) {
        auto currNumBins = possibleNumBins[index];

        // Does this number of bins benefit from this level?
        if (level % ComputeLog(currNumBins) == 0) {
          for (unsigned ptr = 0, limit = histogram[side].first; ptr != limit;
               ++ptr) {
            auto interval = histogram[side].second[ptr];
            // [first, second[ also takes into consideration the `first-1`th
            // element. This is due to `lcp`-array, which takes the previous
            // element into consideration. That's why: `second` - `first` + 1.
            assert(interval.second > interval.first);
            unsigned intervalSize = interval.second - interval.first + 1;

            // Add the length of the interval to all deltas < `intervalSize.
            matrix[index][std::min(intervalSize - 1, maxPossibleTreeError)]
                .first += intervalSize;

            // Does the interval breach the max error, i.e. > max error?
            // Then all deltas < `intervalSize` should receive a `+`.
            if (level != lg) {
              matrix[index][std::min(intervalSize - 1, maxPossibleTreeError)]
                  .second++;
            }
          }
        }
      }
    };

    // Init the histogram.
    ptrInSorted = 0, side = 0, histogram[0].first = 0;
    AddNewInterval({1, spline_points_.size()});

    // And compute the statistics.
    for (unsigned level = 1; level <= lg; ++level) {
      // Init the next row of the histogram.
      side = 1 - side;
      histogram[side].first = 0;

      // And analyze the intervals.
      for (unsigned index = 0, limit = histogram[1 - side].first;
           index != limit; ++index) {
        AnalyzeInterval(level, histogram[1 - side].second[index]);
      }

      // Finally, consume the current level.
      ConsumeLevel(level);
    }

    // Update `statistics`.
    const auto UpdateWith = [&](unsigned index, unsigned delta) -> void {
      statistics.emplace_back(
          possibleNumBins[index], delta,
          1.0 * matrix[index][delta].first / spline_points_.size() +
              ComputeLog(delta, true),
          static_cast<size_t>(1 + matrix[index][delta].second) *
              possibleNumBins[index] * sizeof(unsigned));
    };

    // Build the prefix sums.
    for (unsigned index = 0; index != numPossibleBins; ++index) {
      // Finalize statistics in a backward pass.
      // This pass simply represents the computation of prefix sums.
      UpdateWith(index, maxPossibleTreeError);
      for (unsigned delta = maxPossibleTreeError; delta != 1; --delta) {
        matrix[index][delta - 1].first += matrix[index][delta].first;
        matrix[index][delta - 1].second += matrix[index][delta].second;
        UpdateWith(index, delta - 1);
      }
    }
  }

  void ComputeStatistics(std::vector<Statistics>& statistics) {
    ComputeRadixTableStatistics(statistics);
    ComputeCHTStatistics(statistics);
  }

  Statistics InferTuning(std::vector<Statistics>& statistics) {
    assert(!statistics.empty());

    // Find best cost under the given space limit.
    const size_t space_limit =
        static_cast<size_t>(spline_points_.size()) * sizeof(Coord);
    unsigned bestIndex = 0;
    for (unsigned index = 1, limit = statistics.size(); index != limit;
         ++index) {
      const auto elem = statistics[index];
      if (elem.space > space_limit) continue;
      if ((elem.cost < statistics[bestIndex].cost) ||
          ((std::fabs(elem.cost - statistics[bestIndex].cost) < precision) &&
           (elem.space < statistics[bestIndex].space))) {
        bestIndex = index;
      }
    }
    return statistics[bestIndex];
  }

  const Coord min_;
  const Coord max_;
  const Y spline_max_error_;
  std::vector<Coord> spline_points_;

  size_t curr_num_coords_;
  size_t curr_num_distinct_coords_;
  Coord prev_coord_;
  cht::Builder<X> chtb_;

  // Current upper and lower limits on the error corridor of the spline.
  Coord upper_limit_;
  Coord lower_limit_;
};

}  // namespace invertible_cdf::plex
