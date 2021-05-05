from pycbio.sys.objDict import ObjDict
from pycbio.hgdata.bed import Bed

class Range(ObjDict):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __len__(self):
        return self.end - self.start

class Coords(Range):
    def __init__(self, chrom, start, end, strand):
        super().__init__(start, end)
        self.chrom = chrom
        self.strand = strand

class MappedExon(ObjDict):
    """exon mapped to another assembly.  This will list all of the INDEL in
    source exon. start/end can be None if exon is not mapped"""
    def __init__(self, srcExonId, src, srcBases, mapped=None, mappedBases=0):
        self.srcExonId = srcExonId
        self.src = src
        self.srcBases = srcBases
        self.mapped = mapped
        self.mappedBases = mappedBases

class MappedTranscript(ObjDict):
    """Mapping of a transcripts.

    srcTransId - e.g. ENST00000327381.7
    mappedTransId -  e.g. ENST00000327381.7-1, where -N is used to handle multiple mappings
    """
    def __init__(self, srcTransId, mappedTransId, src, mapped,
                 geneId, geneName, geneType, transcriptName, transcriptType,
                 cds, exons):
        self.srcTransId = srcTransId
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

def mappedTranscriptToBed(trans):
    if trans.mapped is None:
        return None
    blocks = [Bed.Block(e.mapped.start, e.mapped.end)
              for e in trans.exons if e.mapped is not None]
    mt = trans.mapped
    return Bed(mt.chrom, mt.start, mt.end, trans.mappedTransId, score=0, strand=mt.strand,
               thickStart=trans.cds.start, thickEnd=trans.cds.end, itemRgb=0, blocks=blocks)
