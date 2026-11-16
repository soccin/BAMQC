# shellcheck shell=bash
##
## Sourced by the conpair scripts. conpairContainer prints the default
## Conpair container: $CONPAIR_SIF, else the path below. The one place the
## default is set; a script's own container argument or option overrides it.
##

conpairContainer() {
    echo "${CONPAIR_SIF:-/home/soccin/Code/ConPair/conpair_v0.3.3.sif}"
}
