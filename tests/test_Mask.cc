#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include <sstream>

#include "IBDmix/Mask.h"

TEST(Mask, CanTestInMask) {
  std::istringstream mask_input(
      "1 100 120\n"
      "1 130 140\n"
      "1 160 161\n"
      "1 190 200\n"
      "1 260 281\n"
      "2 130 140\n"
      "4 130 140\n"
      "8 130 140\n"
      "chr8 130 140\n");

  Mask mask(&mask_input);

  ASSERT_EQ(mask.get_chromosomes(), (std::vector<std::string>{"1", "2", "4", "8", "chr8"}));

  ASSERT_FALSE(mask.in_mask("1", 90));

  // check same position multiple times
  ASSERT_FALSE(mask.in_mask("1", 90));
  ASSERT_FALSE(mask.in_mask("1", 90));
  ASSERT_FALSE(mask.in_mask("1", 90));

  // boundary: position at start is not in mask (start, end]
  ASSERT_FALSE(mask.in_mask("1", 100));

  // position just after start is in mask
  ASSERT_TRUE(mask.in_mask("1", 101));

  ASSERT_TRUE(mask.in_mask("1", 102));

  // boundary: position at end is in mask (start, end]
  ASSERT_TRUE(mask.in_mask("1", 120));

  // skip a range
  ASSERT_TRUE(mask.in_mask("1", 161));

  ASSERT_TRUE(mask.in_mask("1", 261));

  // skip chromosomes
  ASSERT_FALSE(mask.in_mask("2", 130));

  ASSERT_FALSE(mask.in_mask("3", 130));

  ASSERT_TRUE(mask.in_mask("4", 131));

  ASSERT_FALSE(mask.in_mask("5", 131));

  ASSERT_FALSE(mask.in_mask("6", 131));

  ASSERT_TRUE(mask.in_mask("8", 131));

  ASSERT_TRUE(mask.in_mask("chr8", 131));

  ASSERT_FALSE(mask.in_mask("chr8", 145));

  // chromosome not in mask
  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));

  // can query backward (since we're loading everything into memory)
  ASSERT_TRUE(mask.in_mask("chr8", 131));
}

TEST(Mask, StandardChromosomeOrdered) {
  std::istringstream mask_input(
    "1 20 30\n"
    "1 40 50\n"
    "2 140 150\n"
    "4 240 250\n"
    "5 260 270\n"
  );

  Mask mask(&mask_input);

  ASSERT_EQ(mask.get_chromosomes(), (std::vector<std::string>{"1", "2", "4", "5"}));

  ASSERT_TRUE(mask.in_mask("1", 25));
  ASSERT_TRUE(mask.in_mask("1", 45));
  ASSERT_FALSE(mask.in_mask("3", 100));
  ASSERT_TRUE(mask.in_mask("4", 245));
  ASSERT_TRUE(mask.in_mask("5", 265));
}

TEST(Mask, DifficultChromosomeOrdered) {
  // THIS IS THE KEY TEST: Numeric chromosome ordering where "8" > "10" lexically
  // but 8 < 10 numerically. The old Mask_Reader failed this test for non-seekable streams.
  // The new Mask should handle this correctly!
  std::istringstream mask_input(
      "8 130 140\n"
      "8 150 160\n"
      "10 380 390\n"
      "10 400 410\n"
  );

  Mask mask(&mask_input);
  ASSERT_EQ(mask.get_chromosomes(), (std::vector<std::string>{"8", "10"}));

  // chr7 not in mask
  ASSERT_FALSE(mask.in_mask("7", 100));

  // chr8 regions work
  ASSERT_TRUE(mask.in_mask("8", 135));
  ASSERT_TRUE(mask.in_mask("8", 155));

  // chr9 not in mask
  ASSERT_FALSE(mask.in_mask("9", 200));

  // chr10 regions work - THIS FIXES THE BUG!
  ASSERT_TRUE(mask.in_mask("10", 385));
  ASSERT_TRUE(mask.in_mask("10", 405));
}

TEST(Mask, EmptyMask) {
  Mask mask;  // Default constructor creates empty mask

  ASSERT_FALSE(mask.in_mask("1", 100));
  ASSERT_FALSE(mask.in_mask("chr1", 200));
  ASSERT_EQ(mask.get_chromosomes().size(), 0);
}

TEST(Mask, NullStreamMask) {
  Mask mask(nullptr);  // nullptr stream creates empty mask

  ASSERT_FALSE(mask.in_mask("1", 100));
  ASSERT_FALSE(mask.in_mask("chr1", 200));
  ASSERT_EQ(mask.get_chromosomes().size(), 0);
}

TEST(Mask, StartGreaterThanEndThrowsFirstLine) {
  std::istringstream mask_input(
    "1 150 140\n"  // start > end
  );

  ASSERT_THROW(Mask mask(&mask_input), std::invalid_argument);
}

TEST(Mask, StartGreaterThanEndThrowsSecondLine) {
  std::istringstream mask_input(
    "1 140 150\n"
    "1 160 150\n"  // start > end on second line
  );

  ASSERT_THROW(Mask mask(&mask_input), std::invalid_argument);
}

TEST(Mask, UnsortedIntervalsThrows) {
  std::istringstream mask_input(
    "1 140 150\n"
    "1 100 120\n"  // out of order: 100 < 140
  );

  ASSERT_THROW(Mask mask(&mask_input), std::invalid_argument);
}

TEST(Mask, NonContiguousChromosomeThrows) {
  std::istringstream mask_input(
    "1 100 120\n"
    "2 200 220\n"
    "1 130 140\n"  // chr1 appears again after chr2
  );

  ASSERT_THROW(Mask mask(&mask_input), std::invalid_argument);
}

TEST(Mask, MalformedLineThrows) {
  std::istringstream mask_input(
    "1 100\n"  // missing end position
  );

  ASSERT_THROW(Mask mask(&mask_input), std::invalid_argument);
}

TEST(Mask, OverlappingIntervals) {
  // Test that overlapping intervals work correctly
  std::istringstream mask_input(
    "1 100 200\n"
    "1 150 250\n"  // overlaps with previous
  );

  Mask mask(&mask_input);

  ASSERT_TRUE(mask.in_mask("1", 120));   // in first interval
  ASSERT_TRUE(mask.in_mask("1", 170));   // in both intervals
  ASSERT_TRUE(mask.in_mask("1", 220));   // in second interval
  ASSERT_FALSE(mask.in_mask("1", 260));  // after both
}

TEST(Mask, BoundaryConditions) {
  std::istringstream mask_input(
    "1 100 200\n"
  );

  Mask mask(&mask_input);

  // Interval is (start, end] - start excluded, end included
  ASSERT_FALSE(mask.in_mask("1", 100));  // at start - NOT in mask
  ASSERT_TRUE(mask.in_mask("1", 101));   // just after start - in mask
  ASSERT_TRUE(mask.in_mask("1", 199));   // just before end - in mask
  ASSERT_TRUE(mask.in_mask("1", 200));   // at end - in mask
  ASSERT_FALSE(mask.in_mask("1", 201));  // just after end - NOT in mask
}

TEST(Mask, GetIntervalCount) {
  std::istringstream mask_input(
    "1 100 120\n"
    "1 130 140\n"
    "2 200 220\n"
  );

  Mask mask(&mask_input);

  ASSERT_EQ(mask.get_interval_count("1"), 2);
  ASSERT_EQ(mask.get_interval_count("2"), 1);
  ASSERT_EQ(mask.get_interval_count("3"), 0);  // chromosome not in mask
}

TEST(Mask, RandomAccessQueries) {
  // Test that we can query in any order (backward, forward, random)
  std::istringstream mask_input(
    "1 100 120\n"
    "2 200 220\n"
    "3 300 320\n"
  );

  Mask mask(&mask_input);

  // Query in reverse order
  ASSERT_TRUE(mask.in_mask("3", 310));
  ASSERT_TRUE(mask.in_mask("2", 210));
  ASSERT_TRUE(mask.in_mask("1", 110));

  // Query forward again
  ASSERT_TRUE(mask.in_mask("1", 115));
  ASSERT_TRUE(mask.in_mask("2", 215));
  ASSERT_TRUE(mask.in_mask("3", 315));

  // Random order
  ASSERT_TRUE(mask.in_mask("2", 205));
  ASSERT_TRUE(mask.in_mask("1", 105));
  ASSERT_TRUE(mask.in_mask("3", 305));
}
