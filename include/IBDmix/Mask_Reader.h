#pragma once

#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

class Mask_Reader {
 public:
  explicit Mask_Reader(std::istream *mask) : mask(mask) { readline(); }
  bool in_mask(const std::string &chrom, uint64_t position);

 private:
  std::string chromosome = "";
  std::string prev_chromosome = "";
  uint64_t start, end;
  uint64_t prev_start = 0;
  std::istream *mask = nullptr;
  void readline();
};
