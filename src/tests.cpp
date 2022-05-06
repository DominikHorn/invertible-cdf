#include <gtest/gtest.h>

#include <algorithm>
#include <invertible-cdf.hpp>
#include <limits>
#include <random>

const size_t test_dataset_size = 100'000;
using Key = std::uint64_t;
using KeyLims = std::numeric_limits<Key>;

/**
 * Generates an unsorted, seeded, uniform random dataset over the entire key
 * domain
 *
 * @param dataset_size size of the resulting dataset
 * @param seed rng seed. Defaults to 42
 *
 */
std::vector<Key> generate_unsorted_dataset(size_t dataset_size,
                                           size_t seed = 42) {
  // use merseene twister for pseudo random number generation
  std::mt19937 rng(seed);

  // distribute over entire key domain. Note that this might cause issues with
  // implementations that rely on sentinel values etc. Adress if necessary
  std::uniform_int_distribution<Key> dist(KeyLims::min(), KeyLims::max());

  // reserve + push_back appears to be the fastest way to achieve this
  std::vector<Key> res;
  res.reserve(dataset_size);
  for (size_t i = 0; i < dataset_size; i++) res.push_back(dist(rng));

  return res;
}

/**
 * Generates a sorted, seeded, uniform random dataset over the entire key
 * domain
 *
 * @param dataset_size size of the resulting dataset
 * @param seed rng seed. Defaults to 42
 *
 */
std::vector<Key> generate_sorted_dataset(size_t dataset_size,
                                         size_t seed = 42) {
  // simply rely on unsorted generation internally and sort after the fact
  auto ds = generate_unsorted_dataset(dataset_size, seed);
  std::sort(ds.begin(), ds.end());
  return ds;
}

/// Test whether building on the same data returns always returns the same
/// result, especially across different construction interfaces
TEST(ICDF, BuildReproducible) {
  // obtain a random test dataset
  auto keys = generate_sorted_dataset(test_dataset_size);

  // don't fit() `a` right away
  invertible_cdf::InvertibleCDF<Key> icdf_a;
  invertible_cdf::InvertibleCDF<Key> icdf_b(keys.begin(), keys.end());

  // since only `b` has been fitted to the data they should not match!
  EXPECT_NE(icdf_a, icdf_b);

  // only fit `a` now
  icdf_a.fit(keys.begin(), keys.end());

  // equality after construction
  EXPECT_EQ(icdf_a, icdf_b);

  // idempotency of construction
  icdf_a.fit(keys.begin(), keys.end());
  EXPECT_EQ(icdf_a, icdf_b);
}

/// Test whether equality operator is implemented correctly.
/// PREREQUISITE: `BuildReproducible` must be working!
TEST(ICDF, Equality) {
  // obtain a random test dataset
  auto keys = generate_sorted_dataset(test_dataset_size);

  invertible_cdf::InvertibleCDF<Key> icdf_a, icdf_b, icdf_c;

  icdf_a.fit(keys.begin(), keys.end());
  icdf_b.fit(keys.begin(), keys.end());
  icdf_c.fit(keys.begin(), keys.end());

  // basic equality
  EXPECT_EQ(icdf_a, icdf_b);
  EXPECT_EQ(icdf_b, icdf_c);

  // symmetry
  EXPECT_EQ(icdf_b, icdf_a);
  EXPECT_EQ(icdf_c, icdf_b);

  // transitivity
  EXPECT_EQ(icdf_a, icdf_c);
  EXPECT_EQ(icdf_c, icdf_a);
}

/// Test building on unsorted data
TEST(ICDF, BuildUnsorted) {
  // obtain a random test dataset
  auto keys = generate_unsorted_dataset(test_dataset_size);

  // index unsorted data
  invertible_cdf::InvertibleCDF<Key> unsorted_icdf;
  unsorted_icdf.fit(keys.begin(), keys.end());

  // sort and index on sorted data
  std::sort(keys.begin(), keys.end());
  invertible_cdf::InvertibleCDF<Key> sorted_icdf;
  sorted_icdf.fit(keys.begin(), keys.end());

  // if everything is implemented correctly, unsorted_icdf
  // and sorted_icdf are identical
  EXPECT_EQ(sorted_icdf, unsorted_icdf);
}

/// Test obtaining positions for keys
TEST(ICDF, PosForKey) {
  // obtain a random test dataset
  auto keys = generate_sorted_dataset(test_dataset_size);

  // index data
  invertible_cdf::InvertibleCDF<Key> icdf;
  icdf.fit(keys.begin(), keys.end());

  // invariant: For all keys, their search bound must contain them.
  for (size_t i = 0; i < keys.size(); i++) {
    const auto& key = keys[i];
    const auto bounds = icdf.pos_for_key(key);

    EXPECT_LE(bounds.min, i);
    EXPECT_GE(bounds.max, i);
  }
}

/// Test obtaining keys for positions gives correct bounds
TEST(ICDF, KeyForPosBounds) {
  // obtain a random test dataset
  auto keys = generate_sorted_dataset(test_dataset_size);

  // index data
  invertible_cdf::InvertibleCDF<Key> icdf;
  icdf.fit(keys.begin(), keys.end());

  // invariant: For positions, their keys must map to them.
  for (size_t i = 0; i < keys.size(); i++) {
    const auto& key = keys[i];
    const auto bounds = icdf.key_for_pos(i);

    if (key < bounds.min) {
      std::cout << "min delta: "
                << std::max(bounds.min, key) - std::min(bounds.min, key)
                << std::endl;
    }
    if (key > bounds.max) {
      std::cout << "max delta: "
                << std::max(bounds.max, key) - std::min(bounds.max, key)
                << std::endl;
    }

    EXPECT_LE(bounds.min, key);
    EXPECT_GE(bounds.max, key);
  }
}

/// Test that keys for pos are monotonically increasing
TEST(ICDF, KeyForPosMonotone) {
  // obtain a random test dataset
  auto keys = generate_sorted_dataset(test_dataset_size);

  // index data
  invertible_cdf::InvertibleCDF<Key> icdf;
  icdf.fit(keys.begin(), keys.end());

  // invariant: For positions, their keys must map to them.
  invertible_cdf::Bounds<Key> last_bounds{KeyLims::min(), KeyLims::min()};
  for (size_t i = 0; i < keys.size(); i++) {
    const auto bounds = icdf.key_for_pos(i);

    EXPECT_LE(bounds.min, bounds.max);
    EXPECT_LE(last_bounds.min, bounds.min);
    EXPECT_LE(last_bounds.max, bounds.max);

    last_bounds = bounds;
  }
}
