#include "IBDmix/IBD_Collection.h"

void IBD_Collection::initialize(const Genotype_Reader &reader) {
  IBDs.reserve(reader.get_samples().size());
  for (auto &sample : reader.get_samples())
    IBDs.emplace_back(sample, threshold, &pool, exclusive_end);
}

void IBD_Collection::update(const Genotype_Reader &reader,
                            std::ostream &output) {
  const std::string &chromosome = reader.getChromosome();
  uint64_t position = reader.getPosition();

  // Every base of the genome scores something: its own LOD where a locus exists,
  // lod_prior everywhere else. Charge the bases lying strictly between the
  // previous locus and this one. The gap is the same for every sample, so it is
  // computed once here rather than inside the per-sample loop.
  //
  // Nothing is charged for the first line, across a change of chromosome, or
  // for a repeated or out-of-order position -- multi-allelic rows share a
  // coordinate, and an unsigned subtraction there would wrap into a huge
  // positive score.
  double gap_penalty = 0;
  if (lod_prior != 0 && chromosome == previous_chromosome &&
      position > previous_position + 1) {
    gap_penalty = (position - previous_position - 1) * lod_prior;
  }

  for (unsigned int i = 0; i < IBDs.size(); i++) {
    IBDs[i].add_lod(chromosome, position, reader.getLodScore(i),
                    reader.getLineFilter() | reader.getRecoverType(i), output,
                    gap_penalty);
  }

  previous_chromosome = chromosome;
  previous_position = position;
}

void IBD_Collection::purge(std::ostream &output) {
  for (unsigned int i = 0; i < IBDs.size(); i++) IBDs[i].purge(output);
}

void IBD_Collection::add_recorder(IBD_Collection::Recorder type) {
  switch (type) {
    case IBD_Collection::Recorder::counts:
      for (auto &ibd : IBDs)
        ibd.add_recorder(std::make_shared<CountRecorder>());
      break;
    case IBD_Collection::Recorder::sites:
      for (auto &ibd : IBDs) ibd.add_recorder(std::make_shared<SiteRecorder>());
      break;
    case IBD_Collection::Recorder::lods:
      for (auto &ibd : IBDs) ibd.add_recorder(std::make_shared<LODRecorder>());
      break;
  }
}

void IBD_Collection::writeHeader(std::ostream &strm) const {
  IBDs[0].writeHeader(strm);
}
