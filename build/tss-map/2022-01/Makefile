root = ../..
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
doSrcOrgType: ${srcPsl} ${srcFa}

${srcBed}: ${srcIdCoordFile}
	endSiteFileToBed $< $@.${TMPEXT}
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
destRawPsl = ${mappedDir}/${srcDb}-${destDb}.${type}.raw.pal.gz
destPsl = ${mappedDir}/${srcDb}-${destDb}.${type}.psl.gz
filtStats = ${mappedDir}/${srcDb}-${destDb}.${type}.filtStats
destBed = ${mappedDir}/${srcDb}-${destDb}.${type}.bed

build_pair: ${destBed}

${destRawPsl}: ${srcPsl} ${srcFa} ${synChains}
	@mkdir -p $(dir $@)
	mapGenes ${srcPsl} ${srcFa} ${synChains} ${destDb} $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

${destPsl}: ${destRawPsl} 
	@mkdir -p $(dir $@)
	pslCDnaFilter -maxAligns=1 -statsOut=${filtStats} ${destRawPsl} /dev/stdout \
	    | sort -k 14,14 -k 16,16n |pigz -9c >${destPsl}.${TMPEXT}
	mv -f $@.${TMPEXT} $@

${destBed}: ${destPsl} ${srcBed}
	@mkdir -p $(dir $@)
	endSitePslToBed ${srcBed} $<  $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

endif

clean:
	rm -f data/endsites/*/*.bed data/endsites/*/*.psl data/endsites/*/*.tmp
