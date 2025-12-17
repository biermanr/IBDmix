#include "IBDmix/Mask_Reader.h"

// Initial scan of the mask file to get the chromosome order.
// Detects if the stream is seekable and routes to appropriate scanning method.
void Mask_Reader::detect_seekability_and_scan() {
  if (mask == nullptr) return;

  std::streampos initial_pos = mask->tellg();
  is_seekable = (initial_pos != std::streampos(-1));

  if (is_seekable) {
    perform_full_validation_scan();
  } else {
    perform_minimal_scan_with_warning();
  }
}

// Seekable stream: perform full validation scan and record chromosome order.
// Only performed if the input stream is seekable so that we can reset to the beginning after this 1st pass
// to avoid consuming the mask input.
void Mask_Reader::perform_full_validation_scan() {
  // Seekable stream: perform scan pass and record chromosome order
  while (readline()) {
    if (chromosome == prev_chromosome) continue;

    if (chrom_seen(chromosome)) {
        throw std::invalid_argument("-- IBDmix Error: Mask file not sorted. Chromosome " + chromosome +
                                    " appears multiple times out of order.");
    }

    chromosome_order.push_back(chromosome);
    prev_chromosome = chromosome;
  }

  // reset to beginning
  mask->clear();
  mask->seekg(0);
  chromosome = "";
  prev_chromosome = "";
  readline();
}

// Non-seekable stream: emit warning and read only first line.
// Full validation will be deferred until runtime during in_mask() calls.
void Mask_Reader::perform_minimal_scan_with_warning() {
  // Stream is not seekable (e.g., from process substitution like <(wget ...))
  std::cerr << "-- IBDmix Warning: Mask file stream is not seekable (e.g., from process substitution like <(wget ...)).\n";
  readline();
}


bool Mask_Reader::in_mask(const std::string &geno_chrom, uint64_t geno_position) {
  // test if chrom/position is in mask file
  // true if position is in (start, end]
  if (mask == nullptr) return false;

  for (;;) {
    if (chromosome == "") return false;

    if (geno_chrom == chromosome) {
      if (geno_position <= start) return false;
      else if (end < geno_position) readline();
      else return true;  // start < position <= end
    } 
    else {
      if (is_seekable && !chrom_seen(geno_chrom)) return false;
      if (!is_seekable && chromosome > geno_chrom) return false; // assuming mask and genotype sorted same way
      readline(); 
    }
  }
}

bool Mask_Reader::readline() {
  if (mask == nullptr) {
    chromosome = "";
    return false;
  }

  std::string line;
  if(!std::getline(*mask, line)) {
    chromosome = "";
    return false;
  }

  std::istringstream iss(line);
  if (!(iss >> chromosome >> start >> end)) {
    chromosome = "";
    throw std::invalid_argument("Unable to parse line mask file " + line);
  }

  check_start_end_ordering();
  if(!is_seekable) {
    if (chromosome != prev_chromosome) {
      if (chrom_seen(chromosome)) {
        throw std::invalid_argument("-- IBDmix Error: Mask file not sorted. Chromosome " + chromosome +
                                    " appears multiple times out of order.");
      }
      chromosome_order.push_back(chromosome);
      prev_chromosome = chromosome;
    }
  }
  return true;
}

// Check the current state of the mask reader, raising exceptions for invalid conditions
void Mask_Reader::check_start_end_ordering() {

  if (start > end) {
    throw std::invalid_argument("-- IBDmix Error: Mask file has start > end for " + 
                                chromosome + ":" + std::to_string(start)+"-"+std::to_string(end));
  }
  
  if ((chromosome == prev_chromosome) && (start < prev_start)) {
    throw std::invalid_argument("-- IBDmix Error: Mask file not sorted. " + 
                                prev_chromosome + ":" + std::to_string(prev_start)+"-"+std::to_string(prev_end) +
                                " comes before " +
                                chromosome + ":" + std::to_string(start)+"-"+std::to_string(end));
  }

  prev_start = start;
  prev_end = end;
}

bool Mask_Reader::chrom_seen(const std::string &chrom) {
  // NOTE: This could be optimized with a set if needed
  return std::find(chromosome_order.begin(), chromosome_order.end(), chrom) != chromosome_order.end();
}