#include "IBDmix/Mask_Reader.h"

// Initial scan of the mask file to validate sorting and store chromosome order.
// Only performed if the input stream is seekable so that we can reset to the beginning after this 1st pass.
void Mask_Reader::scan_pass() {
  if (mask == nullptr) return;
  
  // Check if the stream is seekable by trying to get current position
  std::streampos initial_pos = mask->tellg();
  if (initial_pos == std::streampos(-1)) {
    // Stream is not seekable (e.g., from process substitution like <(wget ...))
    // Fall back to single-pass mode with a warning
    std::cerr << "-- IBDmix Warning: Mask file stream is not seekable (e.g., from process substitution like <(wget ...)).\n"
              << "--                 Mask input validation will be skipped. To enable validation, save the mask to a regular file first.\n"
              << "--\n";
    readline();
    return;
  }
  
  readline();

  while (chromosome != "") {

    check_state(); 
    // NOTE most of the validation checks CAN and SHOULD happen on the streaming readline() calls in in_mask().
    // The only thing we need from this scan pass is to build chromosome_order and validate sorting. We
    // don't want to duplicate all the checks here and in in_mask().
    if (chromosome != prev_chromosome) {
      chromosome_order.push_back(chromosome);
    }

    prev_chromosome = chromosome;
    prev_start = start;
    prev_end = end;
    readline();
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
    if (chromosome == "") {
      return false;
    } else if (geno_chrom == chromosome) {
      if (geno_position <= start)
        return false;
      else if (end < geno_position)
        readline();
      else
        return true;  // start < position <= end
    } else {
      // If the chrom is not known to be in the mask, return false and do not advance the mask
      // How does this work with non-seekable masks?
      if (!chrom_seen(geno_chrom)) {
        return false;
      }
      readline();
    }
  }
}

void Mask_Reader::readline() {
  if (mask == nullptr) return;

  std::string line;
  if(!std::getline(*mask, line)) {
    chromosome = "";
    return;
  }

  std::istringstream iss(line);
  if (!(iss >> chromosome >> start >> end)) {
    chromosome = "";
    throw std::invalid_argument("Unable to parse line mask file " + line);
  }
}



// Check the current state of the mask reader, raising exceptions for invalid conditions:
void Mask_Reader::check_state() {

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

  if ((chromosome != prev_chromosome) && (chrom_seen(chromosome))) {
      throw std::invalid_argument("-- IBDmix Error: Mask file not sorted. Chromosome " + chromosome +
                                  " appears multiple times out of order.");
    }
}

bool Mask_Reader::chrom_seen(const std::string &chrom) {
  // NOTE: This could be optimized with a set if needed
  return std::find(chromosome_order.begin(), chromosome_order.end(), chrom) != chromosome_order.end();
}