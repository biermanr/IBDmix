#include "IBDmix/Mask_Reader.h"

void Mask_Reader::scan_pass() {
  if (mask == nullptr) return;
  readline();

  while (chromosome != "") {

    if (start > end) {
      throw std::invalid_argument("Mask file has start > end for " + 
                                  chromosome + ":" + std::to_string(start)+"-"+std::to_string(end));
    }
    
    if (chromosome == prev_chromosome) {
      if (start < prev_start) {
        throw std::invalid_argument("Mask file not sorted. " + 
                                    prev_chromosome + ":" + std::to_string(prev_start)+"-"+std::to_string(prev_end) +
                                    " comes before " +
                                    chromosome + ":" + std::to_string(start)+"-"+std::to_string(end));
      }
    }

    if(chromosome != prev_chromosome){
      if(std::find(chromosome_order.begin(), chromosome_order.end(), chromosome) == chromosome_order.end()) {
        chromosome_order.push_back(chromosome);
      } else {
        throw std::invalid_argument("Mask file not sorted. Chromosome " + chromosome +
                                    " appears multiple times out of order.");
      }
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

bool Mask_Reader::in_mask(const std::string &chrom, uint64_t position) {
  if (mask == nullptr) return false;

  for (;;) {
    if (chromosome == "") {
      return false;
    } else if (chrom == chromosome) {
      if (position <= start)
        return false;
      else if (end < position)
        readline();
      else
        return true;  // start < position <= end
    } else {
      // If the chromosome is not found in the mask, return false and do not advance the mask
      if (find(chromosome_order.begin(), chromosome_order.end(), chrom) == chromosome_order.end()) {
        return false;
      } 
      readline();
    }
  }
}

void Mask_Reader::readline() {
  if (mask == nullptr) return;
  std::string line;
  if (std::getline(*mask, line)) {
    std::istringstream iss(line);
    if (!(iss >> chromosome >> start >> end)) {
      chromosome = "";
      throw std::invalid_argument("Unable to read mask file " + line);
    }
  } else {
    chromosome = "";
  }
}