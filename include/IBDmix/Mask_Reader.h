#pragma once

#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <string>

class Mask_Reader {
 public:
  explicit Mask_Reader(std::istream *mask): mask(mask) { detect_seekability_and_scan(); }
  bool in_mask(const std::string &chrom, uint64_t position);
  const std::vector<std::string> &getChromosomeOrder() const { return chromosome_order; }

 private:
  std::istream *mask = nullptr;
  bool is_seekable = false;

  std::string chromosome = "";
  std::string prev_chromosome = "";
  std::vector<std::string> chromosome_order;
  uint64_t start, end;
  uint64_t prev_start = 0;
  uint64_t prev_end = 0;

  bool readline();
  void detect_seekability_and_scan();
  void perform_full_validation_scan();
  void perform_minimal_scan_with_warning();
  void check_start_end_ordering();
  bool chrom_seen(const std::string &chrom);
};