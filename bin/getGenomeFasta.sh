#!/bin/bash
##
## Print the reference FASTA on IRIS for a genome build reported by
## getGenomeBuildBAM.sh. Exits 1 when there is no FASTA for the build, when
## its .fai or .dict is missing, or when its .dict is not the build (checked
## with getGenomeBuildBAM.sh, so the FASTA always matches the BAM's
## contigs).
##

set -euo pipefail

if [ "$#" != "1" ]; then
    echo "usage: getGenomeFasta.sh BUILD" >&2
    exit 1
fi

SDIR=$(dirname "$(readlink -f "$0")")
BUILD=$1

BIC=/data1/core001/rsrc/genomic/bic/assemblies
IGENOMES=/data1/core001/rsrc/genomic/mskcc-igenomes/igenomes
REFERENCES=/data1/core001/rsrc/genomic/references

case "$BUILD" in
    b37)                      FASTA=$IGENOMES/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta ;;
    hg19)                     FASTA=$BIC/H.sapiens/hg19/hg19.fasta ;;
    hg19-mainOnly)            FASTA=$IGENOMES/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa ;;
    GRCh38)                   FASTA=$IGENOMES/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta ;;
    GRCh38_GDC/GRCh38.d1.vd1) FASTA=$BIC/H.sapiens/GRCh38.d1.vd1/GRCh38.d1.vd1.fa ;;
    hg38)                     FASTA=$BIC/H.sapiens/hg38/hg38.fasta ;;
    hg38-igenomes)            FASTA=$IGENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa ;;
    mm10)                     FASTA=$BIC/M.musculus/mm10/mm10.fasta ;;
    GRCm38)                   FASTA=$REFERENCES/Mus_musculus/BIC_MSK/GRCm38/Sequence/WholeGenomeFasta/GRCm38.fa ;;
    *)
        echo "ERROR: no reference FASTA for build [$BUILD]" >&2
        exit 1
        ;;
esac

FASTA=$(readlink -f "$FASTA")
DICT=${FASTA%.*}.dict
for F in "$FASTA" "$FASTA.fai" "$DICT"; do
    [ -f "$F" ] || { echo "ERROR: not found: $F" >&2; exit 1; }
done

DICT_BUILD=$("$SDIR/getGenomeBuildBAM.sh" "$DICT" || true)
if [ "$DICT_BUILD" != "$BUILD" ]; then
    echo "ERROR: $DICT is build [$DICT_BUILD], not [$BUILD]" >&2
    exit 1
fi

echo "$FASTA"
