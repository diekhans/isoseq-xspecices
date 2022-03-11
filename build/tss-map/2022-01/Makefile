root = ../../..
include ${root}/defs.mk

srcDbs =  hg38 mm39

ifneq ($(filter all%,${type}),)
    dataSet = all
    srcBedType = bed6
    srcBedAsOpt = 
    destBedType = bed12
    destBedAsOpt =
else
    dataSet = expressed
    srcBedType = bed9+2
    srcBedAsOpt = --as=${etcDir}/endsiteSrc.as
    destBedType = bed12+2
    destBedAsOpt = --as=${etcDir}/endsiteMapped.as
endif

# filtered data
endsitesDataDir = data/endsites
filtSrcIdCoordFile = ${endsitesDataDir}/${srcDb}/${type}ID_coord_cellType_count

# pre-filtered data
allSrcDataDir = data/annotatedPeakFiles_CAGE_PolyA
allSrcBed_hg38_alltss = ${allSrcDataDir}/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
allSrcBed_hg38_allpolya = ${allSrcDataDir}/atlas.clusters.2.0.GRCh38.96_chrNames.bed.gz
allSrcBed_mm39_alltss = ${allSrcDataDir}/mm39_CAGEpeaks_liftedOver.bed.gz
allSrcBed_mm39_allpolya = ${allSrcDataDir}/atlas.clusters_chr.mm39_liftedOver_polyA2.bed.gz
allSrcInBed = ${allSrcBed_${srcDb}_${type}}

# common formats
srcDataDir = processed/${srcDb}
srcPrefix = ${srcDataDir}/${srcDb}.${type}
srcBed = ${srcPrefix}.bed.gz
srcLocBed = ${srcPrefix}.loc.bed.gz
srcFa = ${srcPrefix}.fa.gz
srcPsl = ${srcPrefix}.psl.gz
srcBigBed = ${hubSrcDir}/${srcDb}.${type}.bigBed


types = tss polya alltss allpolya
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

ifeq (${dataSet}, all)
${srcBed}: ${allSrcInBed}
	@mkdir -p $(dir $@)
	smtcat $< | tawk 'length($$1) <= 5 {print $$1,$$2,$$3,$$4,0,$$6}'  | pigz -c > $@.${TMPEXTGZ}
	mv -f $@.${TMPEXTGZ} $@

${srcLocBed}: ${srcBed}
	ln -f $< $@

else
${srcBed}: ${srcIdCoordFile}
	@mkdir -p $(dir $@)
	endSiteFileToBed ${srcDb} $< $@.${TMPEXTGZ}
	mv -f $@.${TMPEXTGZ} $@

${srcLocBed}: ${srcBed}
	@mkdir -p $(dir $@)
	smtcat $< | cut -f 1-6 | csort -k 1,1 -k 2,2n -u | pigz -c > $@.${TMPEXTGZ}
	mv -f $@.${TMPEXTGZ} $@
endif


${srcPsl}: ${srcLocBed}
	@mkdir -p $(dir $@)
	bedToPsl ${srcChromSizes} $< stdout | pigz -c > $@.${TMPEXTGZ}
	mv -f $@.${TMPEXTGZ} $@

${srcFa}: ${srcLocBed}
	@mkdir -p $(dir $@)
	twoBitToFa -bed=${srcLocBed} ${srcTwoBit} stdout | pigz -c >$@.${TMPEXTGZ}
	mv -f $@.${TMPEXTGZ} $@

${srcBigBed}: ${srcBed}
	@mkdir -p $(dir $@)
	buildBigBed --bedType=${srcBedType} ${srcBedAsOpt} ${srcDb} $@.${TMPEXT} ${srcBed}
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
destBed = ${mappedDir}/${srcDb}-${destDb}.${type}.bed.gz
destBigBed = ${hubDestDir}/${srcDb}-${destDb}.${type}.bigBed

build_pair: ${destBed} ${destBigBed}

${destPsl}: ${srcPsl} ${srcFa} ${synChains}
	@mkdir -p $(dir $@)
	mapEndSites ${srcPsl} ${srcFa} ${synChains} ${destDb} $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@

ifeq (${dataSet}, all)
${destBed}: ${destPsl}
	@mkdir -p $(dir $@)
	pslToBed $< /dev/stdout | pigz -c > $@.${TMPEXTGZ}
	mv -f $@.${TMPEXTGZ} $@
else
${destBed}: ${destPsl} ${srcBed}
	@mkdir -p $(dir $@)
	endSitePslToBed ${srcBed} $< /dev/stdout | ${bedSortCmd} > $@.${TMPEXT}
	mv -f $@.${TMPEXT} $@
endif

${destBigBed}: ${destBed}
	buildBigBed --bedType=${destBedType} ${destBedAsOpt} ${destDb} $@.${TMPEXT} ${destBed}
	mv -f $@.${TMPEXT} $@


endif

clean:
	rm -rf processed mapped
