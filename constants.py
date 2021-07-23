from enum import Enum

class Variant_Type(Enum):
    INS = 1
    DEL = 2
    INV = 3
    DUP = 4
    TRANS = 5
    dupINVdup = 6
    delINVdel = 7
    delINVdup = 8
    dupINVdel = 9
    delINV = 10
    INVdel = 11
    dDUP_iDEL = 12
    INS_iDEL = 13
    dupINV = 14
    INVdup = 15
    dDUP = 16

class Operations(Enum):
    DUP = "DUP"
    INV = "INV"
    DEL = "DEL"
    DIS = "DIS"

class Constants():

    # for processing from yaml config file
    MAX_LENGTH_ATTR = "max_length"
    MIN_LENGTH_ATTR = "min_length"
    TYPE_ATTR = "type"
    NUM_ATTR = "number"

    # for Structural Variant class
    SV_KEY = {Variant_Type.INS: [("A",), ("A",)],     
            Variant_Type.DEL: [("A",), ()],
            Variant_Type.INV: [("A",), ("a",)],
            Variant_Type.DUP: [("A",), ("A","A")],
            Variant_Type.TRANS: [("A","_","B"), ("B","_","A")],
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
        