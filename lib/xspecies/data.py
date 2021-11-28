from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.bed import Bed
from pycbio.hgdata.frame import Frame
from pycbio.hgdata import dnaOps

class Region(ObjDict):
    def __init__(self, start, end):
        assert start <= end, f"region: {start} > {end}"
        self.start = start
        self.end = end

    def __len__(self):
        return self.end - self.start

    def overlaps(self, rng):
        return (self.start < rng.end) and (self.end > rng.start)

    def contains(self, rng):
        return (self.start >= rng.start) and (self.end <= rng.end)

    def intersect(self, rng):
        return Region(max(self.start, rng.start),
                      min(self.end, rng.end))

class Coords(Region):
    def __init__(self, chrom, start, end, strand):
        super().__init__(start, end)
        self.chrom = chrom
        self.strand = strand

    def intersect(self, rng):
        n = super().intersect(rng)
        n.chrom = self.chrom
        n.strand = self.strand
        return n

class MappedExon(ObjDict):
    """exon mapped to another assembly.  This will list all of the INDEL in
    source exon. start/end can be None if exon is not mapped"""
    def __init__(self, srcExonId, src, srcBases, *, mapped=None, mappedBases=0,
                 cds=None, frame=None):
        self.srcExonId = srcExonId
        self.src = src
        self.srcBases = srcBases
        self.mapped = mapped
        self.mappedBases = mappedBases
        self.cds = cds
        self.frame = frame

    def __str__(self):
        return f"{self.srcExonId} {self.src} => {self.mapped}"

class MappedTranscript(ObjDict):
    """Mapping of a transcripts.

,    srcTransId - e.g. ENST00000327381.7
    mappedTransId -  e.g. ENST00000327381.7-1, where -N is used to handle multiple mappings
    """
    def __init__(self, srcGenome, srcTransId, mappedGenome, mappedTransId, src, mapped,
                 geneId, geneName, geneType, transcriptName, transcriptType,
                 cds, exons):
        self.srcGenome = srcGenome
        self.srcTransId = srcTransId
        self.mappedGenome = mappedGenome
        self.mappedTransId = mappedTransId
        self.src = src
        self.mapped = mapped
        self.geneId = geneId
        self.geneName = geneName
        self.geneType = geneType
        self.transcriptName = transcriptName
        self.transcriptType = transcriptType
        self.cds = cds
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
    return Bed(mt.chrom, mt.start, mt.end, trans.mappedTransId, score=0, strand=mt.strand,
               thickStart=trans.cds.start, thickEnd=trans.cds.end, itemRgb=0, blocks=blocks)

def _getOverlappedExons(mappedGp, exonRegion):
    return [e for e in mappedGp.exons
            if exonRegion.overlaps(e)]

def _getFrame(strand, overExons):
    # first from start or end with frame
    if strand == '-':
        overExons = reversed(overExons)
    for exon in overExons:
        if exon.frame >= 0:
            return Frame(exon.frame)
    raise Exception(f"frame not found: {overExons}")

def getMappedExonFrame(mappedGp, exonRegion):
    """exonRegion is the mapped region of the source exon, return cdsRegion, frame,
     will all being None if no CDS overlap
    """
    mappedCds = Region(mappedGp.cdsStart, mappedGp.cdsEnd)
    if not exonRegion.overlaps(mappedCds):
        return (None, None)
    overExons = _getOverlappedExons(mappedGp, exonRegion)
    exonCds = mappedCds.intersect(exonRegion)
    return exonCds, _getFrame(mappedGp.strand, overExons)

def getGenomeTwoBit(hgdb):
    return "/hive/data/genomes/{asm}/{asm}.2bit".format(asm=hgdb)

class ChromRegionSeq:
    """used to store partial sequence"""
    def __init__(self, coords, seqreader):
        self.coords = coords
        self.seq = seqreader[coords.chrom][coords.start, coords.end]
        if coords.strand == '-':
            self.seq = dnaOps.reverseComplement(self.seq)

    def get(self, rng):
        assert self.coords.contains(rng)
        start = rng.start - self.coords.start
        return self.seq[start: start + len(rng)]
