#!/bin/bash
##
## Print the genome build of a BAM, CRAM, SAM or sequence dictionary
## (.dict), identified by the md5 of its sorted @SQ SN/LN lines. Prints
## "unknown MD5" and exits 1 for a build not in the table; add it there.
##
## The md5 values come from the BIC getGenomeBuildBAM.sh scripts, which
## sorted under en_US.UTF-8, so the locale is fixed here. Under the C
## locale every value is different.
##

set -euo pipefail

if [ "$#" != "1" ]; then
    echo "usage: getGenomeBuildBAM.sh BAM|CRAM|SAM|DICT" >&2
    exit 1
fi

if ! command -v samtools > /dev/null; then
    echo "ERROR: samtools not on PATH" >&2
    exit 1
fi

## SN and LN are picked by tag, so an @SQ line with another field order
## gives the same hash as the standard SN, LN order.
GENOME_MD5=$(samtools view -H "$1" \
    | awk -F'\t' '$1 == "@SQ" {
        sn = ""; ln = ""
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^SN:/) sn = $i
            else if ($i ~ /^LN:/) ln = $i
        }
        print "@SQ\t" sn "\t" ln
    }' \
    | LC_ALL=en_US.UTF-8 sort \
    | md5sum \
    | awk '{print $1}')

case $GENOME_MD5 in
    # Human Genomes - b37/GRCh37 variants
    b879c678e7fd80718cf20d10c6b846e4)
    # b37 gatk /juno/depot/assemblies/H.sapiens/b37/b37.dict
    # Same SN/LN as igenomes GATK GRCh37 human_g1k_v37_decoy.dict
    echo "b37"
    ;;

    117fce86b797081e0af6d69cbd94dcde)
    # b37 version used by DMP pipeline
    echo "b37_dmp"
    ;;

    41f6cc5fbccb9081d64751271f4afa10)
    # new version of b37 used by DMP new set of decoy's and ???
    echo "b37_dmp_v2"
    ;;

    933b376d936c265fc6b44c8bd19fc66d)
    # TCGA BAMs UR:ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz
    # AS:GRCh37-lite (b37-ish)
    echo "GRCh37-lite"
    ;;

    # Human Genomes - hg19 variants
    5b4e380a6b4fc3494cfc66c917d41b37)
    # UCSC hg19 /juno/depot/assemblies/H.sapiens/hg19/hg19.dict
    echo "hg19"
    ;;

    3d72c6961689390556ed2d5a33e66e17)
    # Main chromosomes only (used by cfDNA collaboration)
    # Same SN/LN as igenomes UCSC hg19 genome.dict
    echo "hg19-mainOnly"
    ;;

    # Human Genomes - GRCh38
    688157c599f72b79f8877cf96bf21fd9)
    # GRCh38 GATK Homo_sapiens_assembly38.dict
    echo "GRCh38"
    ;;

    6a02a9445433ae0b01febe3b4f39d5ec)
    # https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
    echo "GRCh38_GDC/GRCh38.d1.vd1"
    ;;

    492b7d5290a5fd8ff50b325caa8fc1b5)
    # UCSC hg38 /data1/core001/rsrc/genomic/bic/assemblies/H.sapiens/hg38/hg38.dict
    echo "hg38"
    ;;

    75dd63fe37a496aa58f867c70ac422cb)
    # UCSC hg38 from igenomes, Homo_sapiens/UCSC/hg38 genome.dict
    echo "hg38-igenomes"
    ;;

    # Mouse Genomes - mm10 variants
    d660fd17a979374182d3ba8b6d76cac0)
    # UCSC mm10 /juno/depot/assemblies/M.musculus/mm10/mm10.dict
    echo "mm10"
    ;;

    34839afd79d8b772037ed7a0e0a4f9c3)
    # UCSC mm10
    echo "mm10_hBRAF_V600E"
    ;;

    # Mouse Genomes - other
    f9cd233a3d5c9540eece434c65f84f1c)
    # mm9 Full
    echo "mm9Full"
    ;;

    0835b244e0adb20253e5fa1c4ee58ec4)
    # mouse_GRCm38
    echo "GRCm38"
    ;;

    # Xenograft Genomes
    7f8c5323ff7e0ff6b5d40efe53eaf050)
    # BIC Xeno-graph genome
    echo "b37+mm10"
    ;;

    # Other Organisms
    fbd0f3da03b422def6bc1a77e1b44e54)
    # GRCz10 D.rerio (zebrafish)
    echo "GRCz10"
    ;;

    bea6557adb311904bf440322ad8500ff)
    # dm3 Drosophila
    echo "dm3"
    ;;

    723e23ee9c8718cef6013286b924125f)
    # WBcel235/ce11 (C. elegans)
    echo "WBcel235"
    ;;

    *)
    echo "unknown $GENOME_MD5"
    exit 1
    ;;
esac
