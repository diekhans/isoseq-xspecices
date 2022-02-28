
# set root=

SHELL = /bin/bash
export SHELLOPTS=pipefail
.SECONDARY:

export PATH := ${root}/../bin:${PATH}
PPID = $(shell echo $${PPID})
TMPEXT = ${HOSTNAME}.${PPID}.tmp


DBS = hg38 mm39 rheMac10 rn6

chainsDir = ${root}/chains/${srcDb}
synChains = ${chainsDir}/${srcDb}-${destDb}.chains.gz

srcChromSizes = /hive/data/genomes/${srcDb}/chrom.sizes
destChromSizes = /hive/data/genomes/${destDb}/chrom.sizes

srcTwoBit = /hive/data/genomes/${srcDb}/${srcDb}.2bit
destTwoBit = /hive/data/genomes/${destDb}/${destDb}.2bit
