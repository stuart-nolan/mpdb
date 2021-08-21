"""
mpdb: Material Property Data Base as python module
periodicTable.py: module to import and use periodicTable data

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import re

def fW(formula):
    reFW = re.compile(r'([A-Z][a-z]*)(\d*)')
    fWDB = {"Ac": "227",
            "Al": "26.98154",
            "Am": "245",
            "Sb": "121.75",
            "Ar": "39.948",
            "As": "74.9216",
            "At": "210",
            "Ba": "137.34",
            "Bk": "247",
            "Be": "9.01218",
            "Bi": "208.9804",
            "B": "10.81",
            "Br": "79.904",
            "Cd": "112.4",
            "Ca": "40.08",
            "Cf": "251",
            "C": "12.011",
            "Ce": "140.12",
            "Cs": "132.9054",
            "Cl": "35.453",
            "Cr": "51.996",
            "Co": "58.9332",
            "Cu": "63.546",
            "Cm": "247",
            "Db": "262",
            "Dy": "162.5",
            "Es": "252",
            "Er": "167.26",
            "Eu": "151.96",
            "Fm": "257",
            "F": "18.998403",
            "Fr": "223",
            "Gd": "157.25",
            "Ga": "69.72",
            "Ge": "72.59",
            "Au": "196.9665",
            "Hf": "178.49",
            "Hs": "265",
            "He": "4.0026",
            "Ho": "164.9304",
            "H": "1.0079",
            "In": "114.82",
            "I": "126.9045",
            "Ir": "192.22",
            "Fe": "55.847",
            "Kr": "83.8",
            "La": "138.9055",
            "Lr": "262",
            "Pb": "207.2",
            "Li": "6.941",
            "Lu": "174.97",
            "Mg": "24.305",
            "Mn": "54.938",
            "Mt": "265",
            "Md": "258",
            "Hg": "200.59",
            "Mo": "95.94",
            "Bh": "262",
            "Nd": "144.24",
            "Ne": "20.179",
            "Np": "237.0482",
            "Ni": "58.7",
            "Nb": "92.9064",
            "N": "14.00674",
            "No": "259",
            "Os": "190.2",
            "O": "15.9994",
            "Pd": "106.4",
            "P": "30.97376",
            "Pt": "195.09",
            "Pu": "244",
            "Po": "209",
            "K": "39.098",
            "Pr": "140.9077",
            "Pm": "145",
            "Pa": "231.0359",
            "Ra": "226.0254",
            "Rn": "222",
            "Re": "186.207",
            "Rh": "102.9055",
            "Rb": "85.4678",
            "Ru": "101.07",
            "Rf": "261",
            "Sm": "150.4",
            "Sc": "44.9559",
            "Sg": "263",
            "Se": "78.96",
            "Si": "28.086",
            "Ag": "107.868",
            "Na": "22.98977",
            "Sr": "87.4678",
            "S": "32.06",
            "Ta": "180.9479",
            "Tc": "97",
            "Te": "127.6",
            "Tb": "158.9254",
            "Tl": "204.37",
            "Th": "232.0381",
            "Tm": "168.9342",
            "Sn": "118.69",
            "Ti": "47.9",
            "W": "183.5",
            "U": "238.029",
            "V": "50.9414",
            "Xi": "131.3",
            "Yb": "173.04",
            "Y": "88.9059",
            "Zn": "65.38",
            "Zr": "91.33"}
    aW = 0
    for (elem,count) in reFW.findall(formula):
        try:
            count = float(count)
        except:
            count = 1

        if elem in fWDB:
            aW += float(fWDB[elem])*count
    return aW

def elemCount(formula):
    reFW = re.compile(r'([A-Z][a-z]*)(\d*)')
    eC = {}
    for (elem,count) in reFW.findall(formula):
        try:
            count = float(count)
        except:
            count = 1

        if elem in eC.keys():
            eC[elem] += count
        else:
            eC.update({elem: count})
            
    return eC
    
if __name__ == "__main__":
    print(fW('Al2O3')) # fW_Al2O3 = 101.96128
    print(elemCount("Al2O3Al2O3"))
