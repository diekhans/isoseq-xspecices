.SECONDARY:

hgDbs = hg38 rheMac10 mm39 rn6
hubDir = ../../hub

hgwdevUrl = https://hgwdev.gi.ucsc.edu/~markd/gencode/isoseq-xspecies/hub/hub.txt
hgwdevUdcDir = ${TMPDIR}/${USER}_xspecies_hgwdev_udc

makeAnnotMapTrackDb = ../bin/makeAnnotMapTrackDb
makeCrgTrackDb = ../bin/makeCrgTrackDb

srcFiles = $(wildcard *.txt *.html */*.txt */*.html)

all: hub

hub: ${hgDbs:%=%.makeDb}

clean:
	rm -f ${hubDir}/*.{txt,html}  ${hubDir}/*/*.{txt,html}


%.makeDb:
	${MAKE} makeDb hgDb=$*

destFiles = ${srcFiles:%=${hubDir}/%}
destTrackDb = ${hubDir}/${hgDb}/trackDb.txt
destTrackAnnotMapDb = ${hubDir}/${hgDb}/annotmap.trackDb.txt
destTrackDbInput = ${destTrackAnnotMapDb}

destTrackCrgDb = ${hubDir}/${hgDb}/crg.trackDb.txt
ifeq (${hgDb}, hg38)
   destTrackDbInput += ${destTrackCrgDb}
endif
ifeq (${hgDb}, mm39)
   destTrackDbInput += ${destTrackCrgDb}
endif



hgwdevUdcDbDir = ${hgwdevUdcDir}_${hgDb}
makeDb: ${destFiles} ${destTrackDb}
	@rm -rf ${hgwdevUdcDbDir}
	hubCheck -udcDir=${hgwdevUdcDbDir} -genome=${hgDb} ${hgwdevUrl}


${hubDir}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@

${destTrackDb}: ${destTrackDbInput} $(wildcard ${hgDb}/*.trackDb.txt)
	@mkdir -p $(dir $@)
	cat $^ > $@.tmp
	mv -f $@.tmp $@

${destTrackAnnotMapDb}: ${makeAnnotMapTrackDb}
	${makeAnnotMapTrackDb} ${hgDb} 1.0 $@.tmp
	mv -f $@.tmp $@

${destTrackCrgDb}: ${makeCrgTrackDb}
	${makeCrgTrackDb} ${hgDb} 10.0 $@.tmp
	mv -f $@.tmp $@
