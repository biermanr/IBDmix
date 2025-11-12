#pragma once

#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <string>

class Mask_Reader {
 public:
  explicit Mask_Reader(std::istream *mask) : mask(mask) { scan_pass(); }
  bool in_mask(const std::string &chrom, uint64_t position);
  const std::vector<std::string> &getChromosomeOrder() const { return chromosome_order; }

 private:
  std::string chromosome = "";
  std::string prev_chromosome = "";
  std::vector<std::string> chromosome_order;
  uint64_t start, end;
  uint64_t prev_start = 0;
  uint64_t prev_end = 0;
  std::istream *mask = nullptr;
  void scan_pass();
  void readline();
};
