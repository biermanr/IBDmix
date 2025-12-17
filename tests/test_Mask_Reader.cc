#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <streambuf>

#include "IBDmix/Mask_Reader.h"

// Non-seekable stream buffer for testing
// Wraps a stringbuf but makes it non-seekable by making seekoff/seekpos fail
class NonSeekableStreamBuf : public std::streambuf {
private:
  std::stringbuf underlying_buf;

public:
  explicit NonSeekableStreamBuf(const std::string& data) : underlying_buf(data) {}

  // Override to make seeking fail (return -1)
  std::streampos seekoff(std::streamoff off, std::ios_base::seekdir dir,
                         std::ios_base::openmode which = std::ios_base::in) override {
    return std::streampos(-1);
  }

  std::streampos seekpos(std::streampos pos,
                         std::ios_base::openmode which = std::ios_base::in) override {
    return std::streampos(-1);
  }

  // Delegate reading to the underlying buffer
  int underflow() override {
    return underlying_buf.sgetc();
  }

  int uflow() override {
    return underlying_buf.sbumpc();
  }

  std::streamsize xsgetn(char* s, std::streamsize n) override {
    return underlying_buf.sgetn(s, n);
  }
};

// Non-seekable stream for testing
class NonSeekableStream : public std::istream {
private:
  NonSeekableStreamBuf buf;

public:
  explicit NonSeekableStream(const std::string& data) : std::istream(&buf), buf(data) {}
};



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
    "1 20 30\n"
    "1 40 50\n"
    "2 140 150\n"
    "4 240 250\n"
    "5 260 270\n"
  );

  // Creating the mask does a scan pass so all chroms are known
  Mask_Reader mask(&mask_input);
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"1", "2", "4", "5"}));

  ASSERT_TRUE(mask.in_mask("1", 25));
  ASSERT_TRUE(mask.in_mask("1", 45));
  ASSERT_TRUE(mask.in_mask("2", 145));

  // chrom missing, but chrom order known after scan so the mask should NOT advance
  ASSERT_FALSE(mask.in_mask("3", 100));

  ASSERT_TRUE(mask.in_mask("4", 245));
  ASSERT_TRUE(mask.in_mask("5", 265));
}

TEST(MaskReader, NonSeekableChromosomeOrdered) {
  NonSeekableStream mask_input(
    "1 20 30\n"
    "1 40 50\n"
    "2 140 150\n"
    "4 240 250\n"
    "5 260 270\n"
  );

  Mask_Reader mask(&mask_input);

  // No scan, but first line read so chrom 1 known
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"1"}));

  mask.in_mask("1", 45); // forces read second line of chrom 1, so still only chrom 1 known
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"1"}));

  mask.in_mask("2", 145); // forces read chrom 2, so now chrom 2 known
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"1", "2"}));

  // checking chr3 forces read of the chr4 line, but then stops since 4 > 3 (assuming sorted)
  ASSERT_FALSE(mask.in_mask("3", 100));
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"1", "2", "4"}));

  ASSERT_TRUE(mask.in_mask("4", 245)); 
  ASSERT_TRUE(mask.in_mask("5", 265)); 
}

TEST(MaskReader, DifficultChromosomeOrdered) {
  // Difficult situation where numerically 8 < 10 but lexically "8" > "10" 
  std::istringstream mask_input(
      "8 130 140\n"
      "8 150 160\n"
      "10 380 390\n"
      "10 400 410\n"
  );

  Mask_Reader mask(&mask_input);
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"8", "10"}));

  ASSERT_FALSE(mask.in_mask("7", 100));
  ASSERT_TRUE(mask.in_mask("8", 135));
  ASSERT_TRUE(mask.in_mask("8", 155));
  ASSERT_FALSE(mask.in_mask("9", 200));
  ASSERT_TRUE(mask.in_mask("10", 385));
  ASSERT_TRUE(mask.in_mask("10", 405));
}

TEST(MaskReader, NonSeekableDifficultChromosomeOrdered) {
  // Difficult situation where numerically 8 < 10 but lexically "8" > "10" 
  NonSeekableStream mask_input(
    "8 130 140\n"
    "8 150 160\n"
    "10 380 390\n"
    "10 400 410\n"
  );

  // Non-seekable stream, so no scan pass, only first line read so only chr8 known
  Mask_Reader mask(&mask_input);
  ASSERT_EQ(mask.getChromosomeOrder(), (std::vector<std::string>{"8"}));

  // chr7 < chr8, so should return false without advancing which is correct
  ASSERT_FALSE(mask.in_mask("7", 100));
  ASSERT_TRUE(mask.in_mask("8", 135));
  ASSERT_TRUE(mask.in_mask("8", 155));

  // checking chr9 forces read of chr10 AND continues since chr10 < chr9 lexically
  ASSERT_FALSE(mask.in_mask("9", 200));

  // this means that these sites which SHOULD be in the mask are already passed.
  // this isn't the ideal behavior, but it's the best we can do without a scan pass
  // or storing the mask in memory or in a temp file for re-reading.
  ASSERT_FALSE(mask.in_mask("10", 385));
  ASSERT_FALSE(mask.in_mask("10", 405));
}

TEST(MaskReader, StartGreaterThanEndThrowsFirstLine) {
  std::istringstream mask_input(
      "1 120 110\n"
      "1 130 140\n"
  );

  // Creating the mask should throw during the scan pass
  ASSERT_THROW(Mask_Reader mask(&mask_input), std::invalid_argument);
}

TEST(MaskReader, StartGreaterThanEndThrowsNotFirstLine) {
  std::istringstream mask_input(
      "1 120 130\n"
      "1 150 140\n"
  );

  // Creating the mask causes a scan pass which will throw on the invalid second line
  ASSERT_THROW(Mask_Reader mask(&mask_input), std::invalid_argument);
}

TEST(MaskReader, ScanningStartGreaterThanEndThrowsFirstLine) {
  NonSeekableStream mask_input(
    "1 140 130\n"
  );

  // Creating the mask should throw since the first line is invalid
  ASSERT_THROW(Mask_Reader mask(&mask_input), std::invalid_argument);
}

TEST(MaskReader, ScanningStartGreaterThanEndThrowsNotFirstLine) {
  NonSeekableStream mask_input(
    "1 120 130\n"
    "1 150 140\n"
  );

  // Creating the mask DOES NOT cause a scan pass since the stream is non-seekable so no exception is thrown here
  Mask_Reader mask(&mask_input);

  // But reading the second line should throw due to invalid start/end
  ASSERT_THROW(mask.in_mask("1", 155), std::invalid_argument);
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