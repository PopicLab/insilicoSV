from enum import Enum

class Variant_Type(Enum):
    INS = "INS"
    DEL = "DEL"
    INV = "INV"
    DUP = "DUP"
    TRA = "TRA"
    dupINVdup = "dupINVdup"
    delINVdel = "delINVdel"
    delINVdup = "delINVdup"
    dupINVdel = "dupINVdel"
    delINV = "delINV"
    INVdel = "INVdel"
    dDUP_iDEL = "dDUP_iDEL"
    INS_iDEL = "INS_iDEL"
    dupINV = "dupINV"
    INVdup = "INVdup"
    dDUP = "dDUP"
    Custom = "Custom"

class Operations(Enum):
    INS = "INS"
    DUP = "DUP"
    INV = "INV"
    DEL = "DEL"
    TRA = "TRA"

class Symbols():
    DIS = "_" # dispersion event              
    PLACEHOLDER = "-"  # placeholder with length 0

# for processing from yaml config file
MAX_LENGTH_ATTR = "max_length"
MIN_LENGTH_ATTR = "min_length"
TYPE_ATTR = "type"
NUM_ATTR = "number"
TRANSFORM_SOURCE_ATTR = "source"
TRANSFORM_TARGET_ATTR = "target"

# for Structural Variant class
SV_KEY = {Variant_Type.INS: [("-",), ("A",)],     
        Variant_Type.DEL: [("A",), ()],
        Variant_Type.INV: [("A",), ("a",)],
        Variant_Type.DUP: [("A",), ("A","A")],
        Variant_Type.TRA: [("A","_","B"), ("B","_","A")],
        Variant_Type.dupINVdup: [("A","B","C"), ("A","c","b","a","C")],
        Variant_Type.delINVdel: [("A","B","C"), ("b",)],
        Variant_Type.delINVdup: [("A","B","C"), ("c","b","C")],
        Variant_Type.dupINVdel: [("A","B","C"), ("A","b","a")],
        Variant_Type.delINV: [("A","B"), ("b",)],
        Variant_Type.INVdel: [("A","B"), ("a",)],
        Variant_Type.dDUP_iDEL: [("A","_","B"), ("A","_","A")],
        Variant_Type.INS_iDEL: [("A","_","B"), ("_","A")],
        Variant_Type.dupINV: [("A","B"), ("A","b","a")],
        Variant_Type.INVdup: [("A","B"), ("b","a","B")],
        Variant_Type.dDUP: [("A","_"), ("A","_","A")]}
        