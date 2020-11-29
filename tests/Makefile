.SECONDARY:

input = input

test: addGtfExonSourcesTests


addGtfExonSourcesTests: addGtfHg38Mm39Test addGtfHg38RheMac10Test


addGtfHg38Mm39Test: mkdirs
	../bin/addGtfExonSources hg38 ${input}/hg38.src.psl ${input}/hg38-mm39.psl ${input}/hg38-mm39.gtf output/$@.gtf
	diff expected/$@.gtf output/$@.gtf

addGtfHg38RheMac10Test: mkdirs
	../bin/addGtfExonSources hg38 ${input}/hg38.src.psl ${input}/hg38-rheMac10.psl ${input}/hg38-rheMac10.gtf output/$@.gtf
	diff expected/$@.gtf output/$@.gtf


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
