from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.psl import PslTbl, PslReader, dropQueryUniq
from pycbio.hgdata.genePred import GenePredTbl
from pycbio.hgdata.bed import Bed
from pycbio.hgdata.frame import Frame
from pycbio.tsv import TsvReader

class Region(ObjDict):
    def __init__(self, start, end):
        assert start <= end, f"region: {start} > {end}"
        self.start = start
        self.end = end

    def __len__(self):
        if self.end <= self.start:
            return 0
        else:
            return self.end - self.start

    def overlaps(self, rng):
        return (self.start < rng.end) and (self.end > rng.start)

    def contains(self, rng):
        return (rng.start >= self.start) and (rng.end <= self.end)

    def intersect(self, rng):
        start = max(self.start, rng.start)
        end = min(self.end, rng.end)
        if end < start:
            end = start
        return Region(start, end)

    def noneIfEmpty(self):
        return self if len(self) > 0 else None

class RCoords(Region):
    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom  # put first
        self.strand = strand
        super().__init__(start, end)

    def intersect(self, rng):
        n = super().intersect(rng)
        n.chrom = self.chrom
        n.strand = self.strand
        return n

    @classmethod
    def fromRegion(cls, chrom, strand, region):
        return RCoords(chrom, region.start, region.end, strand)

class MappingsData:
    def __init__(self, srcPsl, srcGp, srcMeta,
                 mappedPsl, mappedGp):
        assert srcPsl.qSize == mappedPsl.qSize
        assert len(srcPsl.blocks) == len(srcGp.exons)
        assert len(mappedPsl.blocks) == len(mappedGp.exons)
        assert srcPsl.tStrand == '+'
        assert mappedPsl.tStrand == '+'
        self.srcPsl = srcPsl
        self.srcGp = srcGp
        self.srcMeta = srcMeta
        self.mappedPsl = mappedPsl
        self.mappedGp = mappedGp

def loadSrcMeta(metaTsv):
    return {t.transId: t for t in TsvReader(metaTsv)}

def loadMappings(srcPslFile, srcGenePredFile, srcMetaTsv,
                 mappedPslFile, mappedGenePredFile):
    srcPsls = PslTbl(srcPslFile, qNameIdx=True)
    srcGps = GenePredTbl(srcGenePredFile, buildUniqIdx=True)
    srcMetas = loadSrcMeta(srcMetaTsv)
    mappedGps = GenePredTbl(mappedGenePredFile, buildUniqIdx=True)

    mappingData = []
    for mappedPsl in PslReader(mappedPslFile):
        srcName = dropQueryUniq(mappedPsl.qName)
        mappingData.append(MappingsData(srcPsls.qNameMap[srcName][0],
                                        srcGps.names[srcName],
                                        srcMetas[srcName],
                                        mappedPsl,
                                        mappedGps.names[mappedPsl.qName]))
    return mappingData


class MappedExon(ObjDict):
    """Exon mapped to another assembly.  tart/end can be None if exon is not mapped.
    This is saved to the JSON file """
    def __init__(self, srcExonId, srcExonNum, src, srcBases, *, mapped=None, mappedBases=0,
                 srcCds=None, mappedCds=None, frame=None, dnaAlign=None, protAlign=None):
        self.srcExonId = srcExonId
        self.srcExonNum = srcExonNum
        self.src = src
        self.srcBases = srcBases
        self.srcCds = srcCds
        self.mapped = mapped
        self.mappedBases = mappedBases
        self.mappedCds = mappedCds
        self.frame = Frame.fromFrame(frame)
        self.dnaAlign = dnaAlign
        self.protAlign = protAlign

    def __str__(self):
        return f"{self.srcExonId} {self.src} => {self.mapped}"

class MappedTranscript(ObjDict):
    """Mapping of a transcripts.

    srcTransId - e.g. ENST00000327381.7
    mappedTransId -  e.g. ENST00000327381.7-1, where -N is used to handle multiple mappings
    """
    def __init__(self, srcGenome, srcTransId, mappedGenome, mappedTransId, src, mapped,
                 geneId, geneName, geneType, transcriptName, transcriptType,
                 srcCds, mappedCds, exons):
        self.srcGenome = srcGenome
        self.srcTransId = srcTransId
        self.src = src
        self.srcCds = srcCds.noneIfEmpty()
        self.mappedGenome = mappedGenome
        self.mappedTransId = mappedTransId
        self.mapped = mapped
        self.mappedCds = mappedCds.noneIfEmpty()
        self.geneId = geneId
        self.geneName = geneName
        self.geneType = geneType
        self.transcriptName = transcriptName
        self.transcriptType = transcriptType
        self.exons = exons
        self.validate()

    def validate(self):
        """sanity check that there are no overlaping exons."""
        for iExon0 in range(len(self.exons)):
            for iExon1 in range(iExon0 + 1, len(self.exons)):
                me0 = self.exons[iExon0].mapped
                me1 = self.exons[iExon1].mapped
                if ((me0 is not None) and (me1 is not None) and
                    me0.overlaps(me1)):
                    raise Exception(f"{self.srcTransId} {me1} overlaps {me0}")


def mappedTranscriptToBed(trans):
    if trans.mapped is None:
        return None
    blocks = [Bed.Block(e.mapped.start, e.mapped.end)
              for e in trans.exons if e.mapped is not None]
    mt = trans.mapped
    cds = trans.mappedCds if trans.mappedCds is not None else Region(mt.end, mt.end)
    return Bed(mt.chrom, mt.start, mt.end, trans.mappedTransId, score=0, strand=mt.strand,
               thickStart=cds.start, thickEnd=cds.end, itemRgb=0, blocks=blocks)

def getGenomeTwoBit(hgdb):
    return "/hive/data/genomes/{asm}/{asm}.2bit".format(asm=hgdb)
