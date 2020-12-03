.SECONDARY:

hgDbs = hg38 rheMac10 mm39 rn6
hubDir = ../../hub

hgwdevUrl = https://hgwdev.gi.ucsc.edu/~markd/gencode/isoseq-xspecies/hub/hub.txt
hgdevUdcDir = ${TMPDIR}/${USER}_xspecies_hgwdev_udc

srcFiles = $(wildcard *.txt *.html)

all: hub

hub: ${hgDbs:%=%.makeDb}
	@rm -rf ${hgdevUdcDir}
	hubCheck -udcDir=${hgdevUdcDir} ${hgwdevUrl}

%.makeDb:
	${MAKE} makeDb hgDb=$*

destFiles = ${srcFiles:%=${hubDir}/%}
destTrackDb = ${hubDir}/${hgDb}/trackDb.txt
makeDb: ${destFiles} ${destTrackDb}


${hubDir}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@

${destTrackDb}: ../bin/makeTrackDb
	../bin/makeTrackDb ${hgDb} $@.tmp
	mv -f $@.tmp $@
