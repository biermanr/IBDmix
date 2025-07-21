# Usage

## IBDmix Tool

The main `ibdmix` tool identifies IBD segments in genomic data.

### Basic Usage

```bash
./ibdmix -g genotype_file -s sample_file -r mask_file -o output_file
```

### Command Line Options

#### Required Arguments
- `-g, --genotype FILE` - Input genotype file
- `-o, --output FILE` - Output file (use `-` for stdout)

#### Optional Arguments
- `-s, --sample FILE` - Sample file specifying which samples to analyze
- `-r, --mask FILE` - Mask file specifying regions to exclude
- `-d, --LOD-threshold FLOAT` - LOD score threshold (default: 3.0)
- `-m, --minor-allele-count-threshold INT` - Minor allele count threshold (default: 1)
- `-a, --archaic-error FLOAT` - Archaic genotyping error rate (default: 0.01)
- `-e, --modern-error-max FLOAT` - Modern genotyping error rate (default: 0.002)
- `-c, --modern-error-proportion INT` - Modern error proportion (default: 2)

#### Additional Options
- `--more-stats` - Include additional statistics in output
- `--inclusive-end` - Use inclusive end coordinates
- `--write-snps` - Include SNP information in output
- `--write-lods` - Include LOD scores in output
- `-i, --inclusive-end` - Short form for inclusive end
- `-t, --more-stats` - Short form for more stats
- `-w, --write-snps` - Short form for write SNPs

### Examples

#### Basic IBD Detection
```bash
./ibdmix -g data/genotypes.gz -s data/samples.txt -r data/mask.bed -o results.txt
```

#### With Additional Statistics
```bash
./ibdmix -g data/genotypes.gz -s data/samples.txt --more-stats -o results_detailed.txt
```

#### Using Short Arguments
```bash
./ibdmix -g data/genotypes.gz -s data/samples.txt -itw -o results.txt
```

## Generate GT Tool

The `generate_gt` tool prepares genotype files for IBD analysis.

### Basic Usage

```bash
./generate_gt -a archaic_vcf -m modern_vcf -o output_file
```

### Command Line Options

- `-a, --archaic FILE` - Archaic VCF file
- `-m, --modern FILE` - Modern VCF file  
- `-o, --output FILE` - Output genotype file

### Example

```bash
./generate_gt --archaic data/neanderthal.vcf.gz --modern data/modern_humans.vcf.gz --output genotypes.txt
```

## File Formats

### Genotype File Format
The genotype file contains tab-separated values with:
1. Chromosome
2. Position
3. Archaic allele
4. Modern sample genotypes

### Sample File Format
Simple text file with one sample ID per line.

### Mask File Format
BED format file specifying regions to exclude:
```
chr1    1000    2000
chr1    5000    6000
```

## Output Format

The output contains IBD segments with the following columns:
1. Sample ID
2. Chromosome
3. Start position
4. End position
5. LOD score (if --write-lods specified)
6. Additional statistics (if --more-stats specified)
