//
// Created by Artem Novikov on 17.02.2023.
//

#include "IntegerPartition.h"

IntegerPartition::IntegerPartition(uint32_t number, uint32_t terms_count, uint32_t base)
    : number_(number), terms_count_(terms_count), base_(base) {}

IntegerPartition::Iterator IntegerPartition::begin() const {
  return {number_, terms_count_, base_};
}

IntegerPartition::Iterator IntegerPartition::end() const {
  return {0};
}

//////////////////////////////////////////////////////ITERATOR////////////////////////////////////////////

IntegerPartition::Iterator::pointer IntegerPartition::Iterator::operator->() const {
  return &partition_;
}

IntegerPartition::Iterator::reference IntegerPartition::Iterator::operator*() const {
  return partition_;
}

IntegerPartition::Iterator& IntegerPartition::Iterator::operator++() {
  NextPartition();
  return *this;
}

IntegerPartition::Iterator IntegerPartition::Iterator::operator++(int) {
  auto copy = *this;
  NextPartition();
  return copy;
}

bool IntegerPartition::Iterator::operator==(const IntegerPartition::Iterator& other) const {
  return partition_ == other.partition_;
}

bool IntegerPartition::Iterator::operator!=(const IntegerPartition::Iterator& other) const {
  return partition_ != other.partition_;
}

IntegerPartition::Iterator::Iterator(uint32_t terms_count) {
  partition_.resize(terms_count, 0); /// end partition
}

IntegerPartition::Iterator::Iterator(uint32_t number, uint32_t terms_count, uint32_t base) : Iterator(terms_count) {
  /// No possible partitions
  if (number < base * terms_count || terms_count == 0 || number == 0) { return; }

  std::fill(partition_.begin(), partition_.end(), base);
  partition_[0] = number - base * (terms_count - 1);
}

void IntegerPartition::Iterator::NextPartition() {
  uint32_t accumulated = 0;
  for (auto partition_last = std::next(partition_.begin()); partition_last != partition_.end(); ++partition_last) {
    if (partition_.front() >= *partition_last + 2) {
      --partition_.front();
      ++(*partition_last);

      auto new_value = *partition_last;
      std::fill(std::next(partition_.begin()), partition_last, new_value);
      partition_.front() += accumulated - new_value * (std::distance(partition_.begin(), partition_last) - 1);
      return;
    }
    accumulated += *partition_last;
  }

  /// last partition
  partition_.clear();
}