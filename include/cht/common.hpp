#pragma once

#include <cstddef>
#include <cstdint>

namespace invertible_cdf::plex::cht {

struct SearchBound {
  size_t begin;
  size_t end;  // Exclusive.
};

}  // namespace invertible_cdf::plex::cht
