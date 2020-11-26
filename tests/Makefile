.SECONDARY:


test: addGtfExonSourcesTests



addGtfExonSourcesTests: addGtfExonSourcesTest


addGtfExonSourcesTest: mkdirs
	../bin/addGtfExonSources hg38 input/hg38.src.psl input/hg38-mm39.psl input/hg38-mm39.gtf output/$@.gtf




mkdirs:
	@mkdir -p output
