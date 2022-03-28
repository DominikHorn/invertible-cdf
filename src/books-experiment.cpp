#include <cstdint>
#include <invertible-cdf.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include "src/support/datasets.hpp"

using namespace invertible_cdf;  // NOLINT(google-build-using-namespace)

int main() {  // NOLINT(bugprone-exception-escape)
  using Key = std::uint64_t;

  std::cout << "load" << std::endl;
  const std::vector<Key> ds = dataset::load_cached(dataset::ID::BOOKS, 2000000);

  std::cout << "fit" << std::endl;
  InvertibleCDF<Key> icdf(ds.begin(), ds.end());

  std::cout << "test" << std::endl;
  for (size_t i = 0; i < ds.size(); i++) {
    const auto& key = ds[i];

    const auto pos_bounds = icdf.pos_for_key(key);
    if (i < pos_bounds.min || i > pos_bounds.max) {
      throw std::runtime_error(
          "invalid pos_bounds for key " + std::to_string(key) + " at index " +
          std::to_string(i) + ": [" + std::to_string(pos_bounds.min) + ", " +
          std::to_string(pos_bounds.max) + "]");
    }

    const auto key_bounds = icdf.key_for_pos(i);
    if (key < key_bounds.min || key > key_bounds.max) {
      throw std::runtime_error("invalid key_bounds for index " +
                               std::to_string(i) + " expected key " +
                               std::to_string(key) + ": [" +
                               std::to_string(key_bounds.min) + ", " +
                               std::to_string(key_bounds.max) + "]");
    }
  }

  return 0;
}
