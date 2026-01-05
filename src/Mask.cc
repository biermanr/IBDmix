#include "IBDmix/Mask.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

Mask::Mask(std::istream* stream) {
  if (stream != nullptr) {
    load_from_stream(stream);
  }
}

void Mask::load_from_stream(std::istream* stream) {
  std::string line;
  std::string prev_chromosome = "";
  Interval* prev_interval = nullptr;

  while (std::getline(*stream, line)) {
    std::istringstream iss(line);
    std::string chromosome;
    uint64_t start, end;

    if (!(iss >> chromosome >> start >> end)) {
      throw std::invalid_argument("Unable to parse mask file line: " + line);
    }

    Interval interval{start, end};

    // Validate interval
    if (start > end) {
      throw std::invalid_argument("-- IBDmix Error: Mask file has start > end for " +
                                  chromosome + ":" + std::to_string(start) + "-" +
                                  std::to_string(end));
    }

    // Check if chromosome changed
    if (chromosome != prev_chromosome) {
      // Check that chromosome hasn't appeared before (non-contiguous)
      if (!prev_chromosome.empty() && intervals_.find(chromosome) != intervals_.end()) {
        throw std::invalid_argument("-- IBDmix Error: Mask file not sorted. Chromosome " +
                                    chromosome + " appears multiple times out of order.");
      }
      chromosome_order_.push_back(chromosome);
      prev_chromosome = chromosome;
      prev_interval = nullptr;
    }

    // Validate ordering within chromosome
    if (prev_interval != nullptr) {
      validate_interval(chromosome, interval, prev_interval);
    }

    // Add interval to the map
    intervals_[chromosome].push_back(interval);
    prev_interval = &intervals_[chromosome].back();
  }
}

void Mask::validate_interval(const std::string& chrom, const Interval& interval,
                             const Interval* prev_interval) const {
  if (prev_interval != nullptr && interval.start < prev_interval->start) {
    throw std::invalid_argument("-- IBDmix Error: Mask file not sorted. " +
                                chrom + ":" + std::to_string(prev_interval->start) + "-" +
                                std::to_string(prev_interval->end) + " comes before " +
                                chrom + ":" + std::to_string(interval.start) + "-" +
                                std::to_string(interval.end));
  }
}

bool Mask::in_mask(const std::string& chromosome, uint64_t position) const {
  auto chrom_it = intervals_.find(chromosome);
  if (chrom_it == intervals_.end()) {
    return false;  // Chromosome not in mask
  }

  const auto& intervals = chrom_it->second;
  if (intervals.empty()) {
    return false;
  }

  // Binary search for first interval where end >= position
  auto it = std::lower_bound(intervals.begin(), intervals.end(), position,
      [](const Interval& interval, uint64_t pos) {
        return interval.end < pos;
      });

  // Check if position falls in this interval
  if (it != intervals.end() && it->contains(position)) {
    return true;
  }

  return false;
}

std::vector<std::string> Mask::get_chromosomes() const {
  return chromosome_order_;
}

size_t Mask::get_interval_count(const std::string& chromosome) const {
  auto it = intervals_.find(chromosome);
  if (it == intervals_.end()) {
    return 0;
  }
  return it->second.size();
}
