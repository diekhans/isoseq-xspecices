class ObjDict(dict):
    """Dict object where keys are field names.
    This is useful for JSON by doing:
       json.load(fh, object_hook=ObjDict):
    or more efficiently in Python3:
       json.load(fh, object_pairs_hook=ObjDict):
    """

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)


class Range(ObjDict):
    def __init__(self, start, end):
        self.start = start
        self.end = end

class Coords(Range):
    def __init__(self, chrom, start, end, strand):
        super().__init__(start, end)
        self.chrom = chrom
        self.strand = strand

class MappedExon(ObjDict):
    """exon mapped to another assembly.  This will list all of the INDEL in
    source exon. start/end can be None if exon is not mapped"""
    def __init__(self, srcExonId, src, mapped=None, mappedBases=0):
        self.srcExonId = srcExonId
        self.src = src
        self.mapped = mapped
        self.mappedBases = mappedBases

class MappedTranscript(ObjDict):
    """Mapping of a transcripts.

    srcTransId - e.g. ENST00000327381.7
    mappedTransId -  e.g. ENST00000327381.7-1, where -N is used to handle multiple mappings
    """
    def __init__(self, srcTransId, mappedTransId, src, mapped, exons):
        self.srcTransId = srcTransId
        self.mappedTransId = mappedTransId
        self.src = src
        self.mapped = mapped
        self.exons = exons
