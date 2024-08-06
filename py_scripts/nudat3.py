import json
import numpy as np
import pint
import contextlib
import os 

home = os.environ['HOME']

elements_dict = {
    1: "h",
    2: "he",
    3: "li",
    4: "be",
    5: "b",
    6: "c",
    7: "n",
    8: "o",
    9: "f",
    10: "ne",
    11: "na",
    12: "mg",
    13: "al",
    14: "si",
    15: "p",
    16: "s",
    17: "cl",
    18: "ar",
    19: "k",
    20: "ca",
    21: "sc",
    22: "ti",
    23: "v",
    24: "cr",
    25: "mn",
    26: "fe",
    27: "co",
    28: "ni",
    29: "cu",
    30: "zn",
    31: "ga",
    32: "ge",
    33: "as",
    34: "se",
    35: "br",
    36: "kr",
    37: "rb",
    38: "sr",
    39: "y",
    40: "zr",
    41: "nb",
    42: "mo",
    43: "tc",
    44: "ru",
    45: "rh",
    46: "pd",
    47: "ag",
    48: "cd",
    49: "in",
    50: "sn",
    51: "sb",
    52: "te",
    53: "i",
    54: "xe",
    55: "cs",
    56: "ba",
    57: "la",
    58: "ce",
    59: "pr",
    60: "nd",
    61: "pm",
    62: "sm",
    63: "eu",
    64: "gd",
    65: "tb",
    66: "dy",
    67: "ho",
    68: "er",
    69: "tm",
    70: "yb",
    71: "lu",
    72: "hf",
    73: "ta",
    74: "w",
    75: "re",
    76: "os",
    77: "ir",
    78: "pt",
    79: "au",
    80: "hg",
    81: "tl",
    82: "pb",
    83: "bi",
    84: "po",
    85: "at",
    86: "rn",
    87: "fr",
    88: "ra",
    89: "ac",
    90: "th",
    91: "pa",
    92: "u",
    93: "np",
    94: "pu",
    95: "am",
    96: "cm",
    97: "bk",
    98: "cf",
    99: "es",
    100: "fm",
    101: "md",
    102: "no",
    103: "lr",
    104: "rf",
    105: "db",
    106: "sg",
    107: "bh",
    108: "hs",
    109: "mt",
    110: "ds",
    111: "rg",
    112: "cn",
    113: "nh",
    114: "fl",
    115: "mc",
    116: "lv",
    117: "ts",
    118: "og",
}

periodic_table = {
    "h": 1,
    "he": 2,
    "li": 3,
    "be": 4,
    "b": 5,
    "c": 6,
    "n": 7,
    "o": 8,
    "f": 9,
    "ne": 10,
    "na": 11,
    "mg": 12,
    "al": 13,
    "si": 14,
    "p": 15,
    "s": 16,
    "cl": 17,
    "ar": 18,
    "k": 19,
    "ca": 20,
    "sc": 21,
    "ti": 22,
    "v": 23,
    "cr": 24,
    "mn": 25,
    "fe": 26,
    "co": 27,
    "ni": 28,
    "cu": 29,
    "zn": 30,
    "ga": 31,
    "ge": 32,
    "as": 33,
    "se": 34,
    "br": 35,
    "kr": 36,
    "rb": 37,
    "sr": 38,
    "y": 39,
    "zr": 40,
    "nb": 41,
    "mo": 42,
    "tc": 43,
    "ru": 44,
    "rh": 45,
    "pd": 46,
    "ag": 47,
    "cd": 48,
    "in": 49,
    "sn": 50,
    "sb": 51,
    "te": 52,
    "i": 53,
    "xe": 54,
    "cs": 55,
    "ba": 56,
    "la": 57,
    "ce": 58,
    "pr": 59,
    "nd": 60,
    "pm": 61,
    "sm": 62,
    "eu": 63,
    "gd": 64,
    "tb": 65,
    "dy": 66,
    "ho": 67,
    "er": 68,
    "tm": 69,
    "yb": 70,
    "lu": 71,
    "hf": 72,
    "ta": 73,
    "w": 74,
    "re": 75,
    "os": 76,
    "ir": 77,
    "pt": 78,
    "au": 79,
    "hg": 80,
    "tl": 81,
    "pb": 82,
    "bi": 83,
    "po": 84,
    "at": 85,
    "rn": 86,
    "fr": 87,
    "ra": 88,
    "ac": 89,
    "th": 90,
    "pa": 91,
    "u": 92,
    "np": 93,
    "pu": 94,
    "am": 95,
    "cm": 96,
    "bk": 97,
    "cf": 98,
    "es": 99,
    "fm": 100,
    "md": 101,
    "no": 102,
    "lr": 103,
    "rf": 104,
    "db": 105,
    "sg": 106,
    "bh": 107,
    "hs": 108,
    "mt": 109,
    "ds": 110,
    "rg": 111,
    "cn": 112,
    "nh": 113,
    "fl": 114,
    "mc": 115,
    "lv": 116,
    "ts": 117,
    "og": 118,
}

def readJSON(path):
    with open(path, 'r', encoding='utf-8') as file:
        return json.load(file)

def get_decay(nuclide, data, mode):
    try:
        modes = data[nuclide]["levels"][0]["decayModes"]
    except IndexError:
        return False
    except KeyError:
        return False
    if "observed" in modes:
        for entry in modes["observed"]:
            #print(entry)
            if entry["mode"] == mode:
                return entry
    return False

def getT12(nuclide, data):
    level = data[nuclide]["levels"][0]
    if "halflife" in level:
        return level["halflife"]
    return False

def getNucNetData(mode):
    ureg = pint.UnitRegistry()
    PATH = "./decays_half_lives.json"
    DATA = readJSON(PATH)
    for nuclide in DATA:
        #print(nuclide)
        MODE = get_decay(nuclide, DATA, mode)
        if MODE is False:
            continue
        branchingRatio = MODE["value"]
        if branchingRatio < 10:
            continue

        T12 = getT12(nuclide, DATA)
        if T12 is False:
            continue
        halflife = T12["value"]

        halflife_unit = T12["unit"]
        if halflife_unit == "m": # nudat terrible
            halflife_unit = "min"
        elif halflife_unit == "y":
            halflife_unit = "yr"
        elif halflife_unit == "h":
            halflife_unit = "hr"
        elif ureg.parse_expression(halflife_unit).dimensionality == ureg.Quantity(1, 'eV').dimensionality:
            continue

        Z_abbreviation = ""
        A = ""
        for char in nuclide:
            if char.isalpha():
                Z_abbreviation += char
            elif char.isdigit():
                A += char

        halflife_s = (float(halflife) * ureg.parse_expression(halflife_unit)).to(ureg.second)
        #print(f"{nuclide}: Halflife= {halflife} {halflife_unit} SF= {branchingRatio}")
        print("single_rate")
        print("NuDat3 json")
        print("1")
        print(Z_abbreviation.lower() + A)
        print("2")
        print("he4")
        print(elements_dict[periodic_table[Z_abbreviation.lower()] - 2] + str(int(A) - 4))
        print((np.log(2) / (halflife_s.magnitude)) * branchingRatio/100)
        print("")

fw = open(f"{home}/nucnet-tools-code/all_nudat3_alpha_reactions.txt", "w")
with contextlib.redirect_stdout(fw):
    getNucNetData("ɑ")
fw.close()
