foo() {
  fgrep $1 input/hg38-mm39.psl |pslFmt >map.txt
  fgrep $1 input/hg38.src.psl |pslFmt >src.txt
  fgrep $1  output/addGtfExonSourcesTest.gtf  >out.gtf
}

