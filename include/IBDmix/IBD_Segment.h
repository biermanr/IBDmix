#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Segment_Recorders.h"

class IBD_Segment {
 public:
  IBD_Segment(std::string name, double threshold, IBD_Pool *pool,
              bool exclusive_end = true);
  ~IBD_Segment();
  // gap_penalty scores the bases lying strictly between the previous locus and
  // this one; 0 (the default) reproduces stock IBDmix. IBD_Collection computes
  // it once per line, since it is the same for every sample.
  void add_lod(std::string chromosome, uint64_t position, double lod,
               unsigned char bitmask, std::ostream &output,
               double gap_penalty = 0.0);
  void add_recorder(std::shared_ptr<Recorder> recorder);
  void purge(std::ostream &output);
  int size() const;
  void write(std::ostream &strm) const;
  void writeHeader(std::ostream &strm) const;

 private:
  std::string name;
  double threshold;
  IBD_Stack segment;
  IBD_Pool *pool;
  std::vector<std::shared_ptr<Recorder>> recorders;
  std::string chromosome = "";
  bool exclusive_end;

  void add_node(IBD_Node *node, std::ostream &output);
  void initialize_stats();
  void update_stats_recursive(const IBD_Node *node);
  void update_stats(const IBD_Node *node);
  void report_stats(std::ostream &output);
};

std::ostream &operator<<(std::ostream &strm, const IBD_Segment &segment);
