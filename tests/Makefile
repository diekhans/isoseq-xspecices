.SECONDARY:

input = input

test: addGtfExonSourcesTests mappedAlignTests buildAnnotBigPslTests findSmallGapsTest


addGtfExonSourcesTests: addGtfHg38Mm39Test addGtfHg38RheMac10Test

addGtfHg38Mm39Test: mkdirs
	../bin/addGtfExonSources hg38 ${input}/hg38.src.psl ${input}/hg38-mm39.psl ${input}/hg38-mm39.gtf output/$@.gtf
	diff expected/$@.gtf output/$@.gtf

addGtfHg38RheMac10Test: mkdirs
	../bin/addGtfExonSources hg38 ${input}/hg38.src.psl ${input}/hg38-rheMac10.psl ${input}/hg38-rheMac10.gtf output/$@.gtf
	diff expected/$@.gtf output/$@.gtf


mappedAlignTests: mappedAlnHg38Mm39Test mappedAlnHg38RheMac10Test

mappedAlnHg38Mm39Test: mkdirs
	../bin/mappedAlignToAnnot hg38 mm39 ${input}/hg38.src.psl ${input}/hg38.src.tsv ${input}/hg38-mm39.psl ${input}/hg38-mm39.gp --annotBed=output/$@.bed output/$@.json
	diff expected/$@.json output/$@.json
	diff expected/$@.bed output/$@.bed

mappedAlnHg38RheMac10Test: mkdirs
	../bin/mappedAlignToAnnot hg38 mm39 ${input}/hg38.src.psl ${input}/hg38.src.tsv ${input}/hg38-rheMac10.psl ${input}/hg38-rheMac10.gp --annotBed=output/$@.bed output/$@.json
	diff expected/$@.json output/$@.json
	diff expected/$@.bed output/$@.bed

buildAnnotBigPslTests: buildAnnotBigPslTest buildAnnotSrcTest

buildAnnotBigPslTest: mkdirs
	../bin/buildAnnotBigPsl input/hg38.src.tsv input/hg38-mm39.psl input/hg38.src.cds input/hg38.src.fa /hive/data/genomes/mm39/chrom.sizes output/$@.bigPsl
	bigBedToBed output/$@.bigPsl output/$@.bed
	diff expected/$@.bed output/$@.bed

buildAnnotSrcTest: mkdirs
	../bin/buildAnnotBigPsl input/hg38.src.tsv input/hg38.src.psl input/hg38.src.cds input/hg38.src.fa /hive/data/genomes/hg38/chrom.sizes output/$@.bigPsl
	bigBedToBed output/$@.bigPsl output/$@.bed
	diff expected/$@.bed output/$@.bed

findSmallGapsTest: mkdirs
	../bin/findSmallGaps input/hg38-rheMac10.gp output/$@.bed
	diff expected/$@.bed output/$@.bed


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
