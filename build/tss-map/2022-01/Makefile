root = ../../..
include ${root}/defs.mk

srcDbs =  hg38 mm39

endsitesDataDir = data/endsites
srcDataDir = ${endsitesDataDir}/${srcDb}
srcPrefix = ${srcDataDir}/${srcDb}.${type}
srcIdCoordFile = ${srcDataDir}/${type}ID_coord_cellType_count
srcBed = ${srcPrefix}.bed
srcLocBed = ${srcPrefix}.loc.bed
srcFa = ${srcPrefix}.fa
srcPsl = ${srcPrefix}.psl
srcBigBed = ${hubSrcDir}/${srcDb}.${type}.bigBed


types = tss polya
.PHONEY: do%

all: srcData mappings

# ================================================================
srcData: ${srcDbs:%=doSrcOrg_%}

doSrcOrg_%:
	${MAKE} doSrcOrg srcDb=$*

# srcDb=
doSrcOrg: ${types:%=doSrcOrgType_%}

doSrcOrgType_%:
	${MAKE} doSrcOrgType srcDb=${srcDb} type=$*

# srcDb= type=
doSrcOrgType: ${srcPsl} ${srcFa} ${srcBigBed}

${srcBed}: ${srcIdCoordFile}
	endSiteFileToBed $< /dev/stdout | ${bedSortCmd} > $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

${srcLocBed}: ${srcBed}
	cut -f 1-6 $< | sort -k 1,1 -k 2,2n -u > $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

${srcPsl}: ${srcLocBed}
	bedToPsl ${srcChromSizes} $< $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

${srcFa}: ${srcLocBed}
	twoBitToFa -bed=${srcLocBed} ${srcTwoBit} $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

${srcBigBed}: ${srcBed}
	bedToBigBed -type=bed6+2 -as=${etcDir}/endsiteSrc.as ${srcBed} ${srcChromSizes} $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

# ================================================================
mappedDir = mapped/${destDb}

mappings: ${srcDbs:%=%_build_src}
%_build_src: srcData
	${MAKE} build_dest srcDb=$*

build_dest: ${DBS:%=%_build_dest}
%_build_dest:
	${MAKE} build_type srcDb=${srcDb} destDb=$*

build_type: ${types:%=%_build_type}
%_build_type:
	${MAKE} build_pair srcDb=${srcDb} destDb=${destDb} type=$*


ifeq (${srcDb},${destDb})
build_pair:  # does nothing
else
##
# make mappings
##
mappedDir = mapped/${destDb}
destPsl = ${mappedDir}/${srcDb}-${destDb}.${type}.psl.gz
filtStats = ${mappedDir}/${srcDb}-${destDb}.${type}.filtStats
destBed = ${mappedDir}/${srcDb}-${destDb}.${type}.bed
destBigBed = ${hubDestDir}/${srcDb}-${destDb}.${type}.bigBed

build_pair: ${destBed} ${destBigBed}

${destPsl}: ${srcPsl} ${srcFa} ${synChains}
	@mkdir -p $(dir $@)
	mapEndSites ${srcPsl} ${srcFa} ${synChains} ${destDb} $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

${destBed}: ${destPsl} ${srcBed}
	@mkdir -p $(dir $@)
	endSitePslToBed ${srcBed} $< /dev/stdout | ${bedSortCmd} > $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

${destBigBed}: ${destBed}
	bedToBigBed -type=bed12+2 -as=${etcDir}/endsiteMapped.as ${destBed} ${destChromSizes} $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@


endif

clean:
	rm -f data/endsites/*/*.bed data/endsites/*/*.psl data/endsites/*/*.tmp
	rm -rf mapped
