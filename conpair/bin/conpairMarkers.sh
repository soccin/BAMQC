# shellcheck shell=bash
##
## Sourced by the conpair scripts. conpairMarkers BUILD prints the Conpair
## marker files in the container for a build reported by
## bin/getGenomeBuildBAM.sh, as a path without the .bed/.txt extension.
## Returns 1 when Conpair has no markers for the build.
##
## The GRCh37 markers name contigs 1, 2, ... and the GRCh38 markers chr1,
## chr2, ..., so hg19 (chr1 on GRCh37 coordinates) has no markers.
##

conpairMarkers() {
    local dir=/usr/bin/conpair/data/markers
    local name=autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8
    case "$1" in
        b37|b37_dmp|b37_dmp_v2|GRCh37-lite)
            echo "$dir/GRCh37.$name" ;;
        GRCh38|GRCh38_GDC/GRCh38.d1.vd1|hg38|hg38-igenomes)
            echo "$dir/GRCh38.$name.liftover" ;;
        *)
            return 1 ;;
    esac
}
