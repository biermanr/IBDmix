#ifndef MASK_H
#define MASK_H

#include <cstdint>
#include <istream>
#include <map>
#include <string>
#include <vector>

// Mask represents a set of genomic intervals to exclude from analysis.
// The mask is loaded entirely into memory on construction for efficient querying.
// Intervals are stored sorted by chromosome and position for binary search lookup.
class Mask {
public:
  // Construct empty mask (no masking applied)
  Mask() = default;

  // Construct mask from input stream
  // Expected format: chromosome start end (whitespace-separated)
  // Validates that intervals are sorted and properly formatted
  // Throws std::invalid_argument on malformed input
  explicit Mask(std::istream* stream);

  // Query if a genomic position falls within any masked interval
  // Returns true if position is in (start, end] for any interval on this chromosome
  // Returns false if chromosome is not in mask or position not in any interval
  bool in_mask(const std::string& chromosome, uint64_t position) const;

  // Get ordered list of chromosomes in the mask (for testing/validation)
  std::vector<std::string> get_chromosomes() const;

  // Compatibility method for existing code
  const std::vector<std::string> getChromosomeOrder() const {
    return get_chromosomes();
  }

  // Get number of intervals for a chromosome (for testing/debugging)
  size_t get_interval_count(const std::string& chromosome) const;

private:
  struct Interval {
    uint64_t start;
    uint64_t end;

    // Check if position is in (start, end]
    bool contains(uint64_t pos) const {
      return start < pos && pos <= end;
    }
  };

  // Intervals grouped by chromosome, sorted by start position
  std::map<std::string, std::vector<Interval>> intervals_;

  // Track chromosome order as they appear in the file
  std::vector<std::string> chromosome_order_;

  // Load intervals from stream with validation
  void load_from_stream(std::istream* stream);

  // Validate interval ordering and ranges
  void validate_interval(const std::string& chrom, const Interval& interval,
                        const Interval* prev_interval) const;
};

#endif // MASK_H
