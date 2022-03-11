from collections import namedtuple


CellTypeColor = namedtuple("CellTypeColor",
                           ("humanCT", "mouseCT", "rgb", "color"))

cellTypeColors = (
    CellTypeColor("InterneuronsCGE_VIP", "InhibNeuron", "250,0,0", "Red"),
    CellTypeColor("AstrocytesProto", "AstrocytesFibrous", "250,150,0", "Orange"),
    CellTypeColor("InterneuronsMGE_PV", "InhibNeuron", "250,0,0", "Red"),
    CellTypeColor("ExciteDG", "ExcitDG", "0,150,0", "Green"),
    CellTypeColor("InterneuronsCGE_LAMP5", "InhibNeuron", "250,0,0", "Red"),
    CellTypeColor("GranuleNB", "GranuleNB", "0,150,0", "Green"),
    CellTypeColor("MFOLs", "MFOLs", "0,200,200", "Cyan"),
    CellTypeColor("OPCs", "COPs", "0,200,200", "Cyan"),
    CellTypeColor("ExciteCA", "ExcitCA", "0,150,0", "Green"),
    CellTypeColor("Microglia", "Microglia", "0,0,250", "DarkBlue"),
    CellTypeColor("InterneuronsMGE_SST", "InhibNeuron", "250,0,0", "Red"),
    CellTypeColor("VascFibro", "Vasc", "0,0,250", "DarkBlue"),
    CellTypeColor("AstrocytesFibrous", "AstrocytesFibrous", "250,150,0", "Orange"),
    CellTypeColor("MOLs", "MOLs", "0,200,200", "Cyan"),
    CellTypeColor("Macrophages", "Macrophages", "0,0,250", "DarkBlue"),
    CellTypeColor("VascEndo", "Vasc", "0,0,250", "DarkBlue"),
    CellTypeColor("Ependymal", "Ependymal", "250,0,250", "Pink"),
    CellTypeColor(None, "InhCajalRetzius", "250,0,0", "Red"),
    CellTypeColor(None, "Progenitors", "0,250,0", "LimeGreen"),
)


humanCellTypeToColor = {c.humanCT: c for c in cellTypeColors if c.humanCT is not None}
mouseCellTypeToColor = {c.mouseCT: c for c in cellTypeColors if c.mouseCT is not None}
