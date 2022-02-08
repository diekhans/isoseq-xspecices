table kmerBed
"BED of kmers associated with exons  "
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "Name of item"
    uint   score;      "Score from 0-1000"
    char[1] strand;    "+ or -"
    string exonCoords; "coordianates of exon"
    int exon5dist;     "distance previous exon, or huge number if none, or -1 if not available"
    int exon3dist;     "distance to next exon, or huge number if none, or -1 if not available"
    )
