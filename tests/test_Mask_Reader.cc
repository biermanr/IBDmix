#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include <sstream>

#include "IBDmix/Mask_Reader.h"

TEST(MaskReader, CanTestInMask) {
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

  Mask_Reader mask(&mask_input);

  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"1", "2", "4", "8", "chr8"}));

  ASSERT_FALSE(mask.in_mask("1", 90));

  // check same position
  ASSERT_FALSE(mask.in_mask("1", 90));
  ASSERT_FALSE(mask.in_mask("1", 90));
  ASSERT_FALSE(mask.in_mask("1", 90));

  ASSERT_FALSE(mask.in_mask("1", 100));

  ASSERT_TRUE(mask.in_mask("1", 101));

  ASSERT_TRUE(mask.in_mask("1", 102));

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

  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));

  // don't retain the last read values on EOF
  ASSERT_FALSE(mask.in_mask("chr8", 131));

  // setting mask to null should short the region check
  Mask_Reader mask2(nullptr);
  ASSERT_FALSE(mask2.in_mask("1", 161));
}

TEST(MaskReader, StandardChromosomeOrdered) {
  std::istringstream mask_input(
      "1 380 390\n"
      "1 400 410\n"
      "2 130 140\n"
      "2 150 160\n"
  );

  Mask_Reader mask(&mask_input);
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"1", "2"}));
}

TEST(MaskReader, CustomChromosomeOrdered) {
  std::istringstream mask_input(
      "2 130 140\n"
      "2 150 160\n"
      "1 380 390\n"
      "1 400 410\n"
  );

  Mask_Reader mask(&mask_input);
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"2", "1"}));
}

TEST(MaskReader, StartGreaterThanEndThrows) {
  std::istringstream mask_input(
      "1 140 130\n"
  );

  // Creating the mask triggers a scan pass that should throw on start > end
  ASSERT_THROW(Mask_Reader mask(&mask_input), std::invalid_argument);
}

TEST(MaskReader, DisorderedPositionThrows) {
  std::istringstream mask_input(
      "1 130 140\n"
      "1 260 281\n"
      "2 130 140\n"
      "4 130 140\n"
      "4 120 130\n"
  );

  // Creating the mask triggers a scan pass that should throw on disordered positions
  ASSERT_THROW(Mask_Reader mask(&mask_input), std::invalid_argument);
}

TEST(MaskReader, DisorderedRepeatedChromosomeThrows) {
  std::istringstream mask_input(
      "1 130 140\n"
      "1 260 281\n"
      "2 130 140\n"
      "1 380 390\n"
  );

  // Creating the mask triggers a scan pass that should throw on disordered repeated chromosomes
  ASSERT_THROW(Mask_Reader mask(&mask_input), std::invalid_argument);
}