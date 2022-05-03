#pragma once

#include <cstddef>
#include <cstdint>
#include <ostream>

namespace invertible_cdf::plex {

template <class X, class Y>
struct Coord {
  X x;
  Y y;

  friend bool operator==(const Coord<X, Y>& a, const Coord<X, Y>& b) {
    return a.x == b.x && a.y == b.y;
  }

  friend std::ostream& operator<<(std::ostream& os, const Coord<X, Y>& c) {
    os << "(" << c.x << ", " << c.y << ")";
    return os;
  }
};

struct SearchBound {
  size_t begin;
  size_t end;  // Exclusive.
};

// A radix config.
struct RadixConfig {
  unsigned shiftBits;
  unsigned prevPrefix;
  unsigned prevSplineIndex;
  double cost;
};

// Statistics.
struct Statistics {
  Statistics() {}

  Statistics(unsigned numBins, unsigned treeMaxError, double cost, size_t space)
      : numBins(numBins),
        treeMaxError(treeMaxError),
        cost(cost),
        space(space) {}

  unsigned numBins;
  unsigned treeMaxError;
  double cost;
  size_t space;
};

}  // namespace invertible_cdf::plex
