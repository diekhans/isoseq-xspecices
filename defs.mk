
# set root=

SHELL = /bin/bash
export SHELLOPTS=pipefail
.SECONDARY:

export PATH := ${root}/bin:${PATH}
PPID = $(shell echo $${PPID})
TMPEXT = ${HOSTNAME}.${PPID}.tmp
TMPEXTGZ = ${HOSTNAME}.${PPID}.tmp.gz

DBS = hg38 mm39 rheMac10 rn6

etcDir = ${root}/etc

chainsDir = ${root}/build/chains/${srcDb}
synChains = ${chainsDir}/${srcDb}-${destDb}.chains.gz

srcChromSizes = /hive/data/genomes/${srcDb}/chrom.sizes
destChromSizes = /hive/data/genomes/${destDb}/chrom.sizes

srcTwoBit = /hive/data/genomes/${srcDb}/${srcDb}.2bit
destTwoBit = /hive/data/genomes/${destDb}/${destDb}.2bit

hubRootDir =  ${root}/../hub
hubSrcDir =  ${hubRootDir}/${srcDb}
hubDestDir =  ${hubRootDir}/${destDb}
downloadDir = ${hubRootDir}/download

bedSortCmd = csort -k1,1 -k2,2n
