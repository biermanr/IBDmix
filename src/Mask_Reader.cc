#include "IBDmix/Mask_Reader.h"

// Initial scan of the mask file to get the chromosome order.
// Only performed if the input stream is seekable so that we can reset to the beginning after this 1st pass
// to avoid consuming the mask input.
void Mask_Reader::scan_pass() {
  if (mask == nullptr) return;
  
  std::streampos initial_pos = mask->tellg();
  bool seekable = (initial_pos != std::streampos(-1));

  if (!seekable) {
    // Stream is not seekable (e.g., from process substitution like <(wget ...))
    std::cerr << "-- IBDmix Warning: Mask file stream is not seekable (e.g., from process substitution like <(wget ...)).\n";
    readline();
    return;
  }
 
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


bool Mask_Reader::in_mask(const std::string &geno_chrom, uint64_t geno_position) {
  if (mask == nullptr) return false;

  for (;;) {
    if (chromosome == "") return false;

    if (geno_chrom == chromosome) {
      if (geno_position <= start) return false;
      else if (end < geno_position) readline();
      else return true;  // start < position <= end
    } else {
      if (!chrom_seen(geno_chrom)) return false;
      prev_chromosome = chromosome;
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