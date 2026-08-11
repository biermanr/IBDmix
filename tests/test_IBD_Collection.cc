#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "IBDmix/Genotype_Reader.h"
#include "IBDmix/IBD_Collection.h"

// The --lod-prior gap term is computed once per line in IBD_Collection::update,
// since it is the same for every sample. These tests pin the three cases where
// no gap may be charged -- the first line, a change of chromosome, and a
// repeated or out-of-order position -- and the one case where it must be.
//
// Every row below gives the archaic sample n1 the ALT homozygote and m1 the
// same, so m1 scores a healthy positive LOD (~+0.60) at each site while the
// other three moderns do not. That is enough for m1 to build a segment we can
// then try to break with a gap.
class Collection : public ::testing::Test {
 protected:
  static std::string row(const std::string &chrom, uint64_t pos) {
    return chrom + "\t" + std::to_string(pos) + "\tA\tT\t2\t2\t0\t0\t0\n";
  }

  static std::string header() {
    return "chrom\tpos\tref\talt\tn1\tm1\tm2\tm3\tm4\n";
  }

  static std::string run(const std::string &gt, double lod_prior) {
    std::istringstream genotype(gt);
    Genotype_Reader reader(&genotype, nullptr, 0.01, 0.002, 2, 1e-200, 1,
                           lod_prior);
    std::istream dummy(nullptr);
    reader.initialize(dummy);

    IBD_Collection ibds(0, true, lod_prior);
    ibds.initialize(reader);

    std::ostringstream output;
    while (reader.update()) ibds.update(reader, output);
    ibds.purge(output);
    return output.str();
  }

  // IBD_Collection writes no header, so every line is a called segment.
  static int calls_for(const std::string &output, const std::string &sample) {
    int n = 0;
    std::istringstream iss(output);
    std::string line;
    while (std::getline(iss, line))
      if (line.rfind(sample + "\t", 0) == 0) ++n;
    return n;
  }
};

TEST_F(Collection, WideGapOnOneChromosomeSplitsASegment) {
  std::string gt = header();
  for (uint64_t pos = 1; pos <= 5; ++pos) gt += row("1", pos);
  gt += row("1", 1000000);

  // Stock behaviour: the empty megabase costs nothing, so it is one call.
  ASSERT_EQ(calls_for(run(gt, 0), "m1"), 1);

  // At -1e-5 per base the ~1 Mb of unscored genome costs -10, more than the
  // +3.0 accumulated to its left, so the segment closes and a new one opens.
  ASSERT_EQ(calls_for(run(gt, -1e-5), "m1"), 2);
}

TEST_F(Collection, NoGapChargedAcrossAChromosomeChange) {
  std::string gt = header();
  for (uint64_t pos = 1; pos <= 5; ++pos) gt += row("1", pos);
  // Same coordinate jump as above, but on a new chromosome, where the distance
  // is meaningless. Without the chromosome guard this would split as above.
  for (uint64_t pos = 1000000; pos <= 1000002; ++pos) gt += row("2", pos);

  ASSERT_EQ(calls_for(run(gt, -1e-5), "m1"), 1);
}

TEST_F(Collection, NoGapChargedForARepeatedPosition) {
  std::string gt = header();
  // Multi-allelic rows share a coordinate. position - previous - 1 would wrap
  // on unsigned subtraction and charge ~1.8e19 bases, erasing every call.
  gt += row("1", 1);
  gt += row("1", 1);
  for (uint64_t pos = 2; pos <= 4; ++pos) gt += row("1", pos);

  ASSERT_EQ(calls_for(run(gt, -1e-5), "m1"), 1);
}

TEST_F(Collection, AdjacentPositionsAreChargedNothing) {
  std::string gt = header();
  for (uint64_t pos = 1; pos <= 200; ++pos) gt += row("1", pos);

  // Consecutive bases have no bases strictly between them, so even a brutal
  // prior leaves a densely covered run completely alone.
  ASSERT_EQ(run(gt, -1.69), run(gt, 0));
}
