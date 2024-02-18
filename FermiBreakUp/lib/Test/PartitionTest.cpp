//
// Created by Artem Novikov on 08.03.2023.
//

#include <gtest/gtest.h>

#include <numeric>
#include "Utilities/IntegerPartition.h"

size_t factorial(size_t num) {
  size_t fact = 1;
  for (size_t i = 2; i <= num; ++i) {
    fact *= i;
  }
  return fact;
}

size_t Cnk(size_t n, size_t k) {
  return factorial(n) / factorial(k) / factorial(n - k);
}

TEST(PartitionTests, SmallOnePaddingTest) {
  const size_t partition_number = 6;
  const size_t partition_parts = 3;

  size_t partition_counter = 0;
  for (const auto& partition : IntegerPartition(partition_number, partition_parts, 1)) {
    auto copy = partition;
    size_t weight = 1;
    while (std::prev_permutation(copy.begin(), copy.end())) {
      ++weight;
    }

    partition_counter += weight;
    ASSERT_EQ(std::accumulate(partition.begin(), partition.end(), 0u), partition_number);
  }

  ASSERT_EQ(partition_counter, Cnk(partition_number - 1, partition_parts - 1));
}

TEST(PartitionTests, BigOnePaddingTest) {
  const size_t partition_number = 20;
  const size_t partition_parts = 7;

  size_t partition_counter = 0;
  for (const auto& partition : IntegerPartition(partition_number, partition_parts, 1)) {
    auto copy = partition;
    size_t weight = 1;
    while (std::prev_permutation(copy.begin(), copy.end())) {
      ++weight;
    }

    partition_counter += weight;
    ASSERT_EQ(std::accumulate(partition.begin(), partition.end(), 0u), partition_number);
  }

  ASSERT_EQ(partition_counter, Cnk(partition_number - 1, partition_parts - 1));
}

TEST(PartitionTests, SmallZeroPaddingTest) {
  const size_t partition_number = 7;
  const size_t partition_parts = 3;

  size_t partition_counter = 0;
  for (const auto& partition : IntegerPartition(partition_number, partition_parts, 0)) {
    auto copy = partition;
    size_t weight = 1;
    while (std::prev_permutation(copy.begin(), copy.end())) {
      ++weight;
    }

    partition_counter += weight;
    ASSERT_EQ(std::accumulate(partition.begin(), partition.end(), 0u), partition_number);
  }

  ASSERT_EQ(partition_counter, Cnk(partition_number + partition_parts - 1, partition_parts - 1));
}

TEST(PartitionTests, BigZeroPaddingTest) {
  const size_t partition_number = 16;
  const size_t partition_parts = 4;

  size_t partition_counter = 0;
  for (const auto& partition : IntegerPartition(partition_number, partition_parts, 0)) {
    auto copy = partition;
    size_t weight = 1;
    while (std::prev_permutation(copy.begin(), copy.end())) {
      ++weight;
    }

    partition_counter += weight;
    ASSERT_EQ(std::accumulate(partition.begin(), partition.end(), 0u), partition_number);
  }

  ASSERT_EQ(partition_counter, Cnk(partition_number + partition_parts - 1, partition_parts - 1));
}