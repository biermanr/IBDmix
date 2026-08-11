#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "IBDmix/Genotype_Reader.h"
#include "IBDmix/IBD_Segment.h"

class IBD_Collection {
 public:
  // lod_prior is the LOD charged to each base pair of the genome that carries
  // no information of its own. 0 reproduces stock IBDmix exactly.
  explicit IBD_Collection(double threshold, bool exclusive_end = true,
                          double lod_prior = 0.0)
      : threshold(threshold),
        exclusive_end(exclusive_end),
        lod_prior(lod_prior) {}
  void initialize(const Genotype_Reader &reader);
  void update(const Genotype_Reader &reader, std::ostream &output);
  void purge(std::ostream &output);

  enum Recorder { counts, sites, lods };
  void add_recorder(IBD_Collection::Recorder type);
  void writeHeader(std::ostream &strm) const;

 private:
  std::vector<IBD_Segment> IBDs;
  double threshold;
  bool exclusive_end;
  double lod_prior;
  std::string previous_chromosome = "";
  uint64_t previous_position = 0;
  IBD_Pool pool;
};
