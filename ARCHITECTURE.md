# IBDmix Architecture

IBDmix detects archaic introgressed segments in modern human genomes by
identifying regions of Identity-By-Descent (IBD) shared between a modern
individual and an archaic genome (e.g., Neanderthal, Denisovan). It uses a
cumulative LOD (log-odds) score approach, scanning the genome site by site and
reporting stretches where accumulated evidence exceeds a threshold.

## Build System

CMake 3.14+, C++11. Uses FetchContent to pull in
[CLI11](https://github.com/CLIUtils/CLI11) (argument parsing) and
[GoogleTest](https://github.com/google/googletest) (unit tests). The build
produces three executables in `build/src/`:

| Executable     | Source              | Purpose                                      |
|----------------|---------------------|----------------------------------------------|
| `generate_gt`  | `generate_gt.cc`    | Merge archaic + modern VCFs into genotype file |
| `ibdmix`       | `main.cc`           | Detect IBD segments from genotype file        |
| `gt_lods`      | `tabulate_lods.cc`  | Diagnostic: per-site LOD score table          |

## Pipeline Overview

```
archaic.vcf ─┐
              ├─ generate_gt ─> genotype.tsv ─> ibdmix ─> ibd_output.tsv ─> summary.sh
modern.vcf  ─┘                                   ^
                                                  │
                                              mask.bed (optional)
                                              samples.txt (optional)
```

1. **generate_gt** merges two VCFs via sorted merge-join on position
2. **ibdmix** reads the genotype file, computes per-site LOD scores, and
   detects IBD segments using a maximum-subarray algorithm
3. **summary.sh** filters and sorts results by LOD and length thresholds

## Source File Map

```
src/
├── generate_gt.cc         Entry point for VCF merge executable
├── main.cc                Entry point for IBD detection executable
├── tabulate_lods.cc       Entry point for LOD diagnostic executable
├── vcf_file.cc            VCF parsing and genotype extraction
├── Genotype_Reader.cc     Reads genotype file, applies masks/filters, computes LODs
├── Mask_Reader.cc         Streaming BED-file intersection for site exclusion
├── Sample_Mapper.cc       Maps genotype file columns to sample indices
├── lod_calculator.cc      LOD score computation (the statistical model)
├── IBD_Collection.cc      Manages one IBD_Segment per modern sample
├── IBD_Segment.cc         Per-sample segment detection (maximum-subarray algorithm)
├── IBD_Stack.cc           Linked-list stack + memory pool for segment nodes
├── Segment_Recorders.cc   Optional per-segment statistics (counts, SNP positions, LODs)
├── summary.sh             Post-processing filter/sort script
└── CMakeLists.txt         Build targets

include/IBDmix/
├── vcf_file.h
├── Genotype_Reader.h
├── Mask_Reader.h
├── Sample_Mapper.h
├── lod_calculator.h
├── IBD_Collection.h
├── IBD_Segment.h
├── IBD_Stack.h
└── Segment_Recorders.h
```

## Key Classes and Relationships

```
main.cc
  Genotype_Reader                 Reads genotype file, orchestrates per-site processing
    ├── Sample_Mapper             Column index <-> sample name mapping
    ├── Mask_Reader               Streaming BED-file site exclusion
    └── LodCalculator             Precomputes LOD cache per site
  IBD_Collection                  One IBD_Segment per modern sample
    ├── IBD_Pool (shared)         Memory pool for IBD_Nodes (slab allocator)
    └── IBD_Segment[]             Per-sample segment tracking
          ├── IBD_Stack           Linked list of IBD_Nodes
          │     └── IBD_Node*     {lod, cumulative_lod, position, bitmask, next}
          └── Recorder[]          Optional statistics collectors
```

All `IBD_Segment` instances share a single `IBD_Pool` for node allocation.
Nodes are recycled globally rather than per-segment, avoiding repeated heap
allocations in the inner loop.

## Stage 1: Genotype Generation (`generate_gt`)

`VCF_File` wraps a single VCF input stream. On construction it reads the
`#CHROM` header to enumerate sample names and pre-allocates buffers. Each call
to `read_line()` parses one VCF record, collapsing diploid GT calls (`A|B` or
`A/B`) to their sum: `0` (hom-ref), `1` (het), `2` (hom-alt), `9` (missing).

The merge logic performs a sorted join on genomic position:
- **Archaic-only site** (archaic pos < modern pos): emit archaic genotype + all-zero modern genotypes, but only if the archaic carries a non-reference allele
- **Shared site** (positions equal): validate matching REF alleles, emit both genotype blocks
- **Modern-only site** (modern pos < archaic pos): skip (archaic would be `0`, uninformative)

**Output format:** Tab-delimited, one row per informative site:
```
chrom  pos  ref  alt  <archaic_sample>  <modern_sample_1>  <modern_sample_2>  ...
```

## Stage 2: IBD Detection (`ibdmix`)

### Per-Site Processing (Genotype_Reader)

For each line in the genotype file:

1. **Mask check** (`Mask_Reader`): streaming BED intersection flags sites in
   excluded regions
2. **Allele frequency** (`find_frequency`): counts non-missing modern alleles.
   Sites where the minor allele count is at or below `minor_allele_cutoff`
   (default 1) are deselected (LOD set to 0)
3. **LOD cache** (`LodCalculator`): precomputes three LOD values (for modern
   genotype 0, 1, 2) given the archaic genotype and allele frequency
4. **Per-sample LODs**: looks up each modern sample's LOD from the cache
5. **Recovery flags**: sites filtered by mask/MAF but showing extreme
   archaic-modern discordance (archaic=2/modern=0 or vice versa) are flagged
   for optional reporting

### LOD Score Model (LodCalculator)

Each site is scored by comparing two hypotheses:

- **H1 (IBD):** the modern sample inherited one haplotype from the archaic
  lineage, so genotype similarity is expected
- **H0 (no IBD):** the modern genotype is drawn independently from the
  population allele frequency

The LOD score is `log10(P(modern_GT | H1, archaic_GT) / P(modern_GT | H0))`.

**Error model:** Both archaic and modern genotypes may have sequencing errors.
The modern error rate is adaptive:
```
modern_error = min(modern_error_max, alt_frequency * modern_error_proportion)
```
so rare variants receive a smaller error rate.

**By archaic genotype:**
- **Archaic = 0 or 2** (homozygous): LOD scores distinguish concordant vs
  discordant modern genotypes. Concordance yields positive LOD (evidence for
  IBD); discordance yields negative LOD.
- **Archaic = 1** (heterozygous): all three modern genotype LODs are negative.
  Het archaic sites are uninformative for IBD and act as a penalty.
- **Archaic = 9** (missing): LOD = 0 for all modern genotypes (site ignored).

### Segment Detection (IBD_Segment + IBD_Stack)

The algorithm is an adaptation of the maximum-subarray problem for genomic
intervals:

1. **Accumulate:** per-sample LOD scores are pushed onto the stack as
   `IBD_Node` objects, maintaining a running cumulative LOD sum
2. **Track maximum:** the node with the highest cumulative LOD is bookmarked as
   the segment `end`
3. **Finalize:** when the cumulative LOD drops below zero, the segment from
   `start` to `end` is complete. If `end`'s cumulative LOD exceeds the
   threshold (default 3.0), the segment is emitted
4. **Reprocess tail:** nodes after `end` that caused the negative sum are
   recycled and reprocessed, correctly handling adjacent IBD segments

```
LOD
 ^        * end (max cumulative)
 |       / \
 |      /   \
 |     /     \
 |----/-------\------> 0
 |   /         \  <- cumulative goes negative, emit [start, end)
 |  * start     *
```

**Output format:** one line per detected segment:
```
ID  chrom  start  end  max_LOD  [optional recorder columns]
```

By default the end position is exclusive (the next position in the genotype
file after the maximum). The `--inclusive-end` flag changes this to the last
position with increasing LOD.

### Optional Recorders (Segment_Recorders)

Enabled via CLI flags, these collect per-segment statistics during detection:

| Flag            | Recorder        | Output                                              |
|-----------------|-----------------|-----------------------------------------------------|
| `--more-stats`  | CountRecorder   | Site tallies: total, positive/negative LOD, mask/MAF filter counts |
| `--write-snps`  | SiteRecorder    | Comma-separated list of positive-LOD positions      |
| `--write-lods`  | LODRecorder     | Comma-separated list of positive-LOD values          |

### Memory Management (IBD_Pool)

`IBD_Pool` is a slab allocator for `IBD_Node` objects. It starts with 1024
nodes and doubles capacity on each reallocation. Nodes are recycled via
`reclaim_node/segment/stack` methods. This avoids per-site heap allocations
across all samples sharing the pool.

## Stage 3: Filtering (`summary.sh`)

A shell script that filters `ibdmix` output by minimum segment length and
minimum LOD score, labels rows with a population name, and sorts by sample ID
then start position. Supports stdin/stdout piping for integration with
compression tools.

## Testing

- **Unit tests:** GoogleTest-based, fetched via CMake FetchContent. Cover
  individual classes (LOD calculations, mask reading, sample mapping, etc.)
- **Acceptance tests:** end-to-end tests comparing pipeline output against
  known-good reference files
- **CI:** GitHub Actions workflows for both unit and acceptance tests
