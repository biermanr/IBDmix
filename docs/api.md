# API Reference

This section provides detailed information about the IBDmix codebase structure and key components.

## Core Components

### IBD_Segment
Represents an IBD segment with start/end positions and associated statistics.

### IBD_Collection
Manages collections of IBD segments for analysis and output.

### Genotype_Reader
Handles reading and parsing of genotype files.

### Mask_Reader
Processes mask files to exclude specific genomic regions.

### Sample_Mapper
Maps sample identifiers to internal representations.

### LOD_Calculator
Computes LOD scores for IBD segment detection.

## File Structure

```
src/
├── main.cc              # Main IBDmix application
├── generate_gt.cc       # Genotype file generator
├── IBD_Segment.cc       # IBD segment implementation
├── IBD_Collection.cc    # IBD collection management
├── Genotype_Reader.cc   # Genotype file parsing
├── Mask_Reader.cc       # Mask file processing
├── Sample_Mapper.cc     # Sample mapping utilities
├── lod_calculator.cc    # LOD score calculations
└── vcf_file.cc         # VCF file handling

include/IBDmix/
├── IBD_Segment.h
├── IBD_Collection.h
├── Genotype_Reader.h
├── Mask_Reader.h
├── Sample_Mapper.h
├── lod_calculator.h
└── vcf_file.h

tests/
├── test_IBD_Segment.cc
├── test_IBD_Collection.cc
├── test_Genotype_Reader.cc
├── test_Mask_Reader.cc
├── test_Sample_Mapper.cc
└── test_lod_calculator.cc
```

## Building and Testing

The project uses CMake for building and CTest for running unit tests. See the [Installation](installation.md) guide for detailed instructions.
