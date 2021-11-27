#!/bin/bash -ex

src=hg38
#dest=mm39
dest=rheMac10

mkdir -p debug
for f in ${src}-${dest}.gp ${src}-${dest}.psl ; do (zfgrep $1 ../build/annot-map/mapped/${dest}/$f.gz >debug/$f )& done
for f in ${src}.src.gp ${src}.src.psl ; do (zfgrep $1 ../build/annot-map/data/${src}/$f.gz >debug/$f )& done
wait

genePredToGtf -utr -source=transMap file debug/${src}-${dest}.gp debug/${src}-${dest}.gtf 

for f in debug/*.psl ; do
    pslFmt $f >$f.txt
done         
for f in debug/*.gp ; do
    genePredFmt $f >$f.txt
done         
