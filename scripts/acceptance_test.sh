#!/bin/bash

set -euo pipefail

test_type=$1

ibdmix="build/src/ibdmix"
generate_gt="build/src/generate_gt"
# to work for local tests use `ibdmix` and `generate_gt`
# executables in IBDmix/build/src folder following README
if [[ ! -f $ibdmix ]]; then
    ibdmix="../../../build/src/ibdmix"
    generate_gt="../../../build/src/generate_gt"
fi

echo "starting $test_type"

url_base="https://zenodo.org/records/15127123/files"

# Download mask file because it gets read in multiple passes.
# This also avoid multiple wget calls
mask_url="$url_base/cell_data_masks_chr20.bed"
mask_file="cell_data_masks_chr20.bed"
if [[ ! -f $mask_file ]]; then
    echo "Downloading $mask_file..."
    wget -qO - "$mask_url" > "$mask_file"
fi

# Download and unzip the genotype file because it gets read in multiple passes.
genotype_url="$url_base/cell_data_outputs_genotype_altai_1kg_20.gz"
genotype_file="cell_data_outputs_genotype_altai_1kg_20.genotype"
if [[ ! -f $genotype_file ]]; then
    echo "Downloading and unzipping $genotype_file..."
    wget -qO - "$genotype_url" | zcat > "$genotype_file"
fi

read_result() {
    wget -qO - $1 | zcat
}

run_genotype() {
    $generate_gt \
        -a <(wget -qO - $1 | zcat) \
        -m <(wget -qO - $2 | zcat) \
        -o -
}

run_genotype_long() {
    $generate_gt \
        --archaic <(wget -qO - $1 | zcat) \
        --modern <(wget -qO - $2 | zcat) \
        --output -
}

run_ibd_pop() {
    $ibdmix \
        -g $1 \
        -s <(wget -qO - $2) \
        -r $3 \
        -d 3.0 \
        -m 1 \
        -a 0.01 \
        -e 0.002 \
        -c 2 \
        -o >( tail -n +2 )
    # the tail removes the header, which was absent in reference for these
}

run_ibd_pop_stream() {
    $ibdmix \
        -g $1 \
        -s <(wget -qO - $2) \
        -r <(cat $3) \
        -d 3.0 \
        -m 1 \
        -a 0.01 \
        -e 0.002 \
        -c 2 \
        -o >( tail -n +2 )
    # the tail removes the header, which was absent in reference for these
}

run_ibd_pop_long() {
    $ibdmix \
        --genotype $1 \
        --sample <(wget -qO - $2) \
        --mask $3 \
        --LOD-threshold 3.0 \
        --minor-allele-count-threshold 1 \
        --archaic-error 0.01 \
        --modern-error-max 0.002 \
        --modern-error-proportion 2 \
        --output >( tail -n +2 )
    # the tail removes the header, which was absent in reference for these
}

run_ibd_all_mask() {
    $ibdmix \
        -g $1 \
        -r $2 \
        --output >( tail -n +2 )
    # the tail removes the header, which was absent in reference for these
}

run_ibd_all_no_mask() {
    $ibdmix \
        -g $1 \
        --output >( tail -n +2 )
    # the tail removes the header, which was absent in reference for these
}

run_ibd_extra() {
    $ibdmix \
        -g $1 \
        -s <(wget -qO - "$url_base/cell_data_samples_GWD.txt") \
        $2 \
        --output >( cat )
}

run_ibd_extra_mask() {
    $ibdmix \
        -g $1 \
        -s <(wget -qO - "$url_base/cell_data_samples_GWD.txt") \
        $2 \
        -r $3 \
        --output >( cat )
}

if [[ $test_type == "gen_20" ]]; then
    resultfile="$url_base/cell_data_outputs_genotype_altai_1kg_20.gz"
    mod="$url_base/cell_data_mod_chr20.vcf.gz"
    arch="$url_base/cell_data_altai_chr20.vcf.gz"

    cmp \
        <(read_result $resultfile) \
        <(run_genotype $arch $mod)

elif [[ $test_type == "populations" ]]; then
    for pop in ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD \
        IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI ; do
        resultfile="$url_base/cell_data_outputs_ibd_raw_altai_1kg_${pop}_20.gz"

        sample="$url_base/cell_data_samples_${pop}.txt"

        # Test with streaming of mask file
        echo "Running ${pop} with stream mask"
        cmp \
            <(read_result "$resultfile") \
            <(run_ibd_pop_stream $genotype_file $sample $mask_file)

        # Test with mask file
        echo "Running ${pop} with mask file"
        cmp \
            <(read_result "$resultfile") \
            <(run_ibd_pop $genotype_file $sample $mask_file)
    done
    # one more for the long args
    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_pop_long $genotype_file $sample $mask_file)

elif [[ $test_type == "all_mask" ]]; then
    resultfile="$url_base/cell_data_outputs_ibd_raw_all_with_mask.gz"

    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_all_mask $genotype_file $mask_file)

elif [[ $test_type == "all_no_mask" ]]; then
    resultfile="$url_base/cell_data_outputs_ibd_raw_all_no_mask.gz"

    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_all_no_mask $genotype_file)

elif [[ $test_type == "extra" ]]; then

    echo "with tab"
    resultfile="$url_base/terminal_tab_genotype.out.gz"
    mod="$url_base/terminal_tab_mod_chr22.vcf.gz"
    arch="$url_base/terminal_tab_AltNea_n10000.vcf.gz"
    cmp \
        <(read_result $resultfile) \
        <(run_genotype_long $arch $mod)

    echo "more stats"
    resultfile="$url_base/cell_data_outputs_ibd_raw_GWD_20_stats.gz"
    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_extra $genotype_file "--more-stats")

    echo "inclusive end"
    resultfile="$url_base/cell_data_outputs_ibd_raw_GWD_20_inclusive.gz"
    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_extra $genotype_file "--inclusive-end")

    echo "with snps"
    resultfile="$url_base/cell_data_outputs_ibd_raw_GWD_20_snps.gz"
    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_extra $genotype_file "--write-snps")

    echo "with lods"
    resultfile="$url_base/cell_data_outputs_ibd_raw_GWD_20_lods.gz"
    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_extra $genotype_file "--write-snps --write-lods")

    echo "short args"
    resultfile="$url_base/cell_data_outputs_ibd_raw_GWD_20_itw.gz"
    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_extra $genotype_file "-itw")

    echo "short args mask"
    resultfile="$url_base/cell_data_outputs_ibd_raw_GWD_20_itw_mask.gz"
    cmp \
        <(read_result "$resultfile") \
        <(run_ibd_extra_mask $genotype_file "-itw" $mask_file)

    echo "short args mask string chroms"
    resultfile="$url_base/cell_data_outputs_ibd_raw_GWD_20_itw_mask.gz"

    str_chrom_mask_file="cell_data_masks_chr20_string.bed"
    awk 'BEGIN {OFS="\t"} {$1 = "chr" $1 ; print $0}' $mask_file > $str_chrom_mask_file


    cmp \
        <(read_result "$resultfile" | awk 'BEGIN {OFS="\t"} NR>1{$2 = "chr" $2} {print $0}' | head) \
        <($ibdmix \
            -g <(awk  'BEGIN {OFS="\t"} NR>1{$1 = "chr" $1} {print $0}' $genotype_file) \
            -s <(wget -qO - "$url_base/cell_data_samples_GWD.txt") \
            -itw \
            -r $str_chrom_mask_file \
            --output >( head ) \
         )

else
    echo "Unknown test type $test_type"
    exit 1
fi

echo "passed $test_type"
