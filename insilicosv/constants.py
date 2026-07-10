from enum import Enum
from dataclasses import dataclass


class VariantType(Enum):
    INS = "INS"

    DEL = "DEL"
    INV = "INV"
    DUP = "DUP"
    mCNV = "mCNV"
    INV_DUP = "INV_DUP"
    DUP_INV = "DUP_INV"

    dDUP = "dDUP"
    INV_dDUP = "INV_dDUP"
    dDUP_INV = "dDUP_INV"
    INV_rTRA = "INV_rTRA"
    nrTRA = "nrTRA"
    rTRA = "rTRA"
    INV_nrTRA = "INV_nrTRA"

    delINV = "delINV"
    INVdel = "INVdel"
    dupINV = "dupINV"
    INVdup = "INVdup"

    INS_iDEL = "INS_iDEL"
    dDUP_iDEL = "dDUP_iDEL"

    dupINVdup = "dupINVdup"
    delINVdel = "delINVdel"
    delINVdup = "delINVdup"
    dupINVdel = "dupINVdel"

    SNP = "SNP"
    INDEL = "INDEL"
    DIVERGENCE = "DIVERGENCE"

    CUSTOM = "Custom"

    trEXP = "trEXP"
    trCON = "trCON"


class Syntax:
    DISPERSION = '_'

    DIVERGENCE = '*'
    MULTIPLE_COPIES = '+'

    ANCHOR_START = '('
    ANCHOR_END = ')'


TR = [VariantType.trCON, VariantType.trEXP]

SV_KEY = {
    VariantType.INS: ((), ("A",)),

    VariantType.DEL: (("A",), ()),
    VariantType.INV: (("A",), ("a",)),
    VariantType.DUP: (("A",), ("A", "A+")),
    VariantType.mCNV: (("A",), ("A+",)),
    VariantType.INV_DUP: (("A",), ("A", "a+")),
    VariantType.DUP_INV: (("A",), ("a", "a+")),

    VariantType.dDUP: (("A", "_"), ("A", "_", "A+")),
    VariantType.INV_dDUP: (("A", "_"), ("A", "_", "a+")),
    VariantType.dDUP_INV: (("A", "_"), ("a", "_", "a+")),
    VariantType.INV_nrTRA: (("A", "_"), ("_", "a")),
    VariantType.nrTRA: (("A", "_"), ("_", "A")),
    VariantType.rTRA: (("A", "_", "B"), ("B", "_", "A")),
    VariantType.INV_rTRA: (("A", "_", "B"), ("b", "_", "a")),

    VariantType.delINV: (("A", "B"), ("b",)),
    VariantType.INVdel: (("A", "B"), ("a",)),
    VariantType.dupINV: (("A", "B"), ("A", "b", "a")),
    VariantType.INVdup: (("A", "B"), ("b", "a", "B")),

    VariantType.INS_iDEL: (("A", "_", "B"), ("_", "A")),
    VariantType.dDUP_iDEL: (("A", "_", "B"), ("A", "_", "A")),

    VariantType.dupINVdup: (("A", "B", "C"), ("A", "c", "b", "a", "C")),
    VariantType.delINVdel: (("A", "B", "C"), ("b",)),
    VariantType.delINVdup: (("A", "B", "C"), ("c", "b", "C")),
    VariantType.dupINVdel: (("A", "B", "C"), ("A", "b", "a")),

    VariantType.SNP: (("A",), ("A*",)),
    VariantType.INDEL: ((), ()),

    VariantType.trCON: ((), ()),
    VariantType.trEXP: ((), ()),
}


@dataclass(order=True, frozen=True)
class Symbol:
    """A symbol denoting a region.

    Examples:

    A
    _1

    """

    name: str

    def __str__(self):
        return self.name
