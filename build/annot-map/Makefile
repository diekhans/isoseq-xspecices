
SHELL = /bin/bash
export SHELLOPTS=pipefail
.SECONDARY:

export PATH := ../../bin:${PATH}
PPID = $(shell echo $${PPID})
TMPEXT = ${HOSTNAME}.${PPID}.tmp

DBS = hg38 mm39 rheMac10 rn6
hubRootDir = ../hub

hg38_gtf = data/hg38/gencode.v35.annotation.gtf.gz
mm39_gtf = data/mm39/Mus_musculus.GRCm39.103.gtf.gz
rheMac10_gtf = data/rheMac10/Macaca_mulatta.Mmul_10.101.chr.gtf.gz
rn6_gtf = data/rn6/Rattus_norvegicus.Rnor_6.0.101.gtf.gz

build: build_src


build_src: ${DBS:%=%_build_src}
%_build_src:
	${MAKE} build_dest srcDb=$*

build_dest: ${DBS:%=%_build_dest}
%_build_dest:
	${MAKE} build_pair srcDb=${srcDb} destDb=$*

ifeq (${srcDb},${destDb})
build_pair:  # does nothing
else
##
# src data
##
srcDir = data/${srcDb}
srcGp = ${srcDir}/${srcDb}.src.gp.gz
srcMeta = ${srcDir}/${srcDb}.src.tsv.gz
srcGtf = ${${srcDb}_gtf}
srcFa = ${srcDir}/${srcDb}.src.fa.gz
srcPsl = ${srcDir}/${srcDb}.src.psl.gz
srcCds = ${srcDir}/${srcDb}.src.cds.gz

${srcMeta}: ${srcGp}
${srcGp}: ${srcGff}
	@mkdir -p $(dir $@)
	convertGtf ${srcGtf} ${srcGp} ${srcMeta}

${srcFa}: ${srcGp}
	@mkdir -p $(dir $@)
	getRnaPred -genePredExt -cdsUpper ${srcDb} ${srcGp} all /dev/stdout | pigz -c > ${srcFa}.${TMPEXT}
	mv -f ${srcFa}.${TMPEXT} ${srcFa}

${srcCds}: ${srcPsl}
${srcPsl}: ${srcGp}
	@mkdir -p $(dir $@)
	genePredToFakePsl ${srcDb} ${srcGp} /dev/stdout >(pigz -c >${srcCds}) | pigz -c > ${srcPsl}.${TMPEXT}
	mv -f ${srcPsl}.${TMPEXT} ${srcPsl}

##
# mapping chains
##
chainsDir = chains/${srcDb}
synChains = ${chainsDir}/${srcDb}-${destDb}.chains.gz

${synChains}:
	@mkdir -p $(dir $@)
	buildMappingChains ${srcDb} ${destDb} $@

##
# make mappings
##
mappedDir = mapped/${destDb}
destRawPsl = ${mappedDir}/${srcDb}-${destDb}.raw.pal.gz
destPsl = ${mappedDir}/${srcDb}-${destDb}.psl.gz
filtStats = ${mappedDir}/${srcDb}-${destDb}.filtStats
destGp = ${mappedDir}/${srcDb}-${destDb}.gp.gz
destGtf = ${mappedDir}/${srcDb}-${destDb}.gtf.gz
destGtfDownload = ${hubRootDir}/downloads/annot/${srcDb}-${destDb}.gtf.gz

${destRawPsl}: ${srcPsl} ${srcFa} ${synChains}
	@mkdir -p $(dir $@)
	mapGenes ${srcPsl} ${srcFa} ${synChains} ${destDb} $@

${destPsl}: ${destRawPsl} 
	@mkdir -p $(dir $@)
	pslCDnaFilter -minQSize=20 -repsAsMatch -ignoreNs -bestOverlap -localNearBest=0.01 -minId=0.35 -minCover=0.15 -statsOut=${filtStats} ${destRawPsl} /dev/stdout \
	    | sort -k 14,14 -k 16,16n |pigz -9c >${destPsl}.${TMPEXT}
	mv -f ${destPsl}.${TMPEXT} ${destPsl}

${destGp}: ${destPsl}
	@mkdir -p $(dir $@)
	transMapPslToGenePred -noBlockMerge ${srcGp} ${destPsl} /dev/stdout | pigz -c > ${destGp}.${TMPEXT}
	mv -f ${destGp}.${TMPEXT} ${destGp}

${destGtf}: ${destGp}
	@mkdir -p $(dir $@)
	genePredToGtf -utr -source=transMap file ${destGp} /dev/stdout \
	    | addGtfExonSources ${srcDb} ${srcPsl} ${destPsl} /dev/stdin /dev/stdout \
	    | pigz -c >${destGtf}.${TMPEXT}
	mv -f ${destGtf}.${TMPEXT} ${destGtf}

${destGtfDownload}: ${destGtf}
	@mkdir -p $(dir $@)
	ln -f $< $@

##
# hub files
##
destChromSizes = /hive/data/genomes/${destDb}/chrom.sizes
hubBigPsl = ${hubRootDir}/${destDb}/${srcDb}-${destDb}.bigPsl

${hubBigPsl}: ${srcMeta} ${destPsl} ${srcCds} ${srcFa} ${destChromSizes}
	@mkdir -p $(dir $@)
	bigPslBuild ${srcMeta} ${destPsl} ${srcCds} ${srcFa} ${destChromSizes} ${hubBigPsl}


##
# main target
##
build_pair: ${srcGp} ${srcMeta} ${srcPsl} ${srcCds} ${srcFa} ${synChains} \
	${destRawPsl} ${destPsl} ${destGp} ${destGtf} ${hubBigPsl} ${destGtfDownload}

endif
