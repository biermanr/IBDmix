#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <math.h>

#include "IBDmix/lod_calculator.h"

double original_cal_lod(int source_gt, int target_gt, double pb, double aerr,
                        double merr, double minesp) {
  double pa = 1 - pb;
  if (source_gt == 0 && target_gt == 0)
    return log10(((1 - aerr) * (1 - merr) + aerr * merr) / pa /
                 (1 - aerr * (1 - aerr)));
  if (source_gt == 0 && target_gt == 1) {
    double temp = 0.5 * ((1 - aerr) * merr + aerr * (1 - merr)) / pb /
                      (1 - aerr * (1 - aerr)) +
                  0.5 * ((1 - aerr) * (1 - merr) + aerr * merr) / pa /
                      (1 - aerr * (1 - aerr));
    return log10(temp);
  }
  if (source_gt == 0 && target_gt == 2)
    return log10(((1 - aerr) * merr + aerr * (1 - merr) + minesp) / pb /
                 (1 - aerr * (1 - aerr) + minesp));
  if (source_gt == 1 && target_gt == 0)
    return -log10(pa * (1 + 2 * aerr * (1 - aerr)));
  if (source_gt == 1 && target_gt == 1)
    return -log10(2 * pa * pb * (1 + 2 * aerr * (1 - aerr)));
  if (source_gt == 1 && target_gt == 2)
    return -log10(pb * (1 + 2 * aerr * (1 - aerr)));
  if (source_gt == 2 && target_gt == 0)
    return log10(((1 - aerr) * merr + aerr * (1 - merr) + minesp) / pa /
                 (1 - aerr * (1 - aerr) + minesp));
  if (source_gt == 2 && target_gt == 1)
    return log10(0.5 * ((1 - aerr) * (1 - merr) + aerr * merr) / pb /
                     (1 - aerr * (1 - aerr)) +
                 0.5 * ((1 - aerr) * merr + aerr * (1 - merr)) / pa /
                     (1 - aerr * (1 - aerr)));
  if (source_gt == 2 && target_gt == 2)
    return log10(((1 - aerr) * (1 - merr) + aerr * merr) / pb /
                 (1 - aerr * (1 - aerr)));
  return 0;
}

TEST(LodCalculator, CanCalculateLod) {
  double pbs[] = {0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.99, 1};
  double aerrs[] = {0.01, 0.02, 0.03};
  double minesps[] = {1e-200, 1e-190, 1e-100};
  char gts[] = {'0', '1', '2', '9'};
  for (double aerr : aerrs)
    for (double minesp : minesps) {
      LodCalculator calc(aerr, 0.002, 2, minesp);
      for (double pb : pbs)
        for (char arch : gts) {
          calc.update_lod_cache(arch, pb, true);
          for (char mod : gts) {
            double read = calc.calculate_lod(mod);
            double orig = original_cal_lod(arch - '0', mod - '0', pb, aerr,
                                           calc.get_modern_error(pb), minesp);
            if (isinf(read) && isinf(orig)) continue;
            ASSERT_TRUE(read == orig || abs((read - orig) / orig) < 0.001);
          }
          calc.update_lod_cache(arch, pb, false);  // not selected
          for (char mod : gts) {
            double read = calc.calculate_lod(mod);
            double orig = original_cal_lod(arch - '0', mod - '0', pb, aerr,
                                           calc.get_modern_error(pb), minesp);
            if ((mod == '0' && arch == '2') || (mod == '2' && arch == '0')) {
              if (isinf(read) && isinf(orig)) continue;
              ASSERT_DOUBLE_EQ(read, orig);
            } else {
              ASSERT_EQ(0, read);
            }
          }
        }
    }
}

TEST(LodCalculator, CanGetModernError) {
  LodCalculator calc(0.01, 0.1, 2);

  ASSERT_DOUBLE_EQ(0.1, calc.get_modern_error(0.4));
  ASSERT_DOUBLE_EQ(0.1, calc.get_modern_error(0.6));
  ASSERT_DOUBLE_EQ(0.1, calc.get_modern_error(0.1));
  ASSERT_DOUBLE_EQ(0.1, calc.get_modern_error(0.9));
  ASSERT_DOUBLE_EQ(0.02, calc.get_modern_error(0.01));
  ASSERT_TRUE(abs((calc.get_modern_error(0.99) - 0.02) / 0.02) < 0.001);
}

// --- --lod-prior fork ------------------------------------------------------
// Stock IBDmix scores an uninformative cell 0, an implicit prior of even odds
// on introgression. The fork replaces that 0 with --lod-prior. The invariant is
// narrow on purpose: cells that stock scores 0 take the prior, and NOTHING else
// may move.

TEST(LodCalculator, LodPriorReplacesExactlyTheUninformativeZeros) {
  double pbs[] = {0.01, 0.05, 0.1, 0.5, 0.9, 0.99};
  char gts[] = {'0', '1', '2', '9'};
  const double prior = -1.69;

  for (double pb : pbs)
    for (char arch : gts)
      for (bool selected : {true, false}) {
        LodCalculator stock(0.01, 0.002, 2, 1e-200);
        LodCalculator fork(0.01, 0.002, 2, 1e-200, prior);
        stock.update_lod_cache(arch, pb, selected);
        fork.update_lod_cache(arch, pb, selected);

        for (char mod : gts) {
          double s = stock.calculate_lod(mod);
          double f = fork.calculate_lod(mod);
          std::string where = std::string("archaic ") + arch + " modern " +
                              mod + " selected " + (selected ? "1" : "0") +
                              " pb " + std::to_string(pb);
          if (s == 0)
            ASSERT_EQ(f, prior) << where;
          else
            ASSERT_EQ(f, s) << where;
        }
      }
}

// The RECOVER_0_2 / RECOVER_2_0 diagnostic mismatches are scored even on a line
// that fails the mask or MAF filter -- they are real, strongly negative evidence
// against IBD, and they are what breaks a run-on segment. The prior must not
// overwrite them. Nothing at --lod-prior 0 can catch this being wrong, since
// there the correct and the clobbered code agree.
TEST(LodCalculator, LodPriorLeavesRecoveryCellsAlone) {
  const double prior = -1.69;
  LodCalculator stock(0.01, 0.002, 2, 1e-200);
  LodCalculator fork(0.01, 0.002, 2, 1e-200, prior);

  // archaic 0 with modern 2, on an unselected line
  stock.update_lod_cache('0', 0.2, false);
  fork.update_lod_cache('0', 0.2, false);
  ASSERT_EQ(fork.calculate_lod('2'), stock.calculate_lod('2'));
  ASSERT_LT(fork.calculate_lod('2'), 0);
  ASSERT_EQ(fork.calculate_lod('0'), prior);
  ASSERT_EQ(fork.calculate_lod('1'), prior);

  // archaic 2 with modern 0, the mirror image
  stock.update_lod_cache('2', 0.2, false);
  fork.update_lod_cache('2', 0.2, false);
  ASSERT_EQ(fork.calculate_lod('0'), stock.calculate_lod('0'));
  ASSERT_LT(fork.calculate_lod('0'), 0);
  ASSERT_EQ(fork.calculate_lod('1'), prior);
  ASSERT_EQ(fork.calculate_lod('2'), prior);
}
