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


// Actual testing starts now
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

TEST(MaskReader, ScanningStartGreaterThanEndThrows) {
  NonSeekableStream mask_input(
    "1 140 130\n"
  );

  // With a non-seekable stream, scan_pass() should print a warning instead of throwing
  // The warning indicates that input validation is skipped for non-seekable streams
  testing::internal::CaptureStderr();
  Mask_Reader mask(&mask_input);
  std::string output = testing::internal::GetCapturedStderr();

  ASSERT_TRUE(output.find("not seekable") != std::string::npos);
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