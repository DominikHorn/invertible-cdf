#pragma once

#include <cstddef>
#include <cstdint>

namespace invertible_cdf::_internal {

// A CDF coordinate.
template <class KeyType>
struct Coord {
  KeyType x;
  double y;

  friend bool operator==(const Coord<KeyType>& a, const Coord<KeyType>& b) {
    return a.x == b.x && a.y == b.y;
  }
};

struct SearchBound {
  size_t begin;
  size_t end;  // Exclusive.
};

}  // namespace invertible_cdf::_internal
