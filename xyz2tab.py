#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#xyz2tab

import sys
import re
#import itertools                                                #for r-length tuples, in sorted order, no repeated elements
#import pandas as pd                                             #pandas tables
#import numpy as np                                              #for calculations
#from scipy.spatial.distance import pdist, squareform, cosine    #for the calculations of the distance matrix and angles (cosine)
#from tabulate import tabulate                                   #nice table output
#import matplotlib.pyplot as plt                                 #for molecule display
                     #for molecule display
#from mpl_toolkits.mplot3d import proj3d                         #for fancy arrows in xyz
from os import path
from pathlib import Path

from lookup_tables import avg_bond_lengths
from helpers import parse_args
from print_tables import PrintTab
from create_plots import CreatePlot

#pd.set_option("display.max_rows", None, "display.max_columns", None)

def classify_bond(pair, bond_length):
    [atom1, atom2] = pair.split("–")

    match = re.search(r"^([a-zA-Z]+)(\d+)", atom1)
    atom1_label = match.group(1)
    atom1_id = match.group(2)
    match = re.search(r"^([a-zA-Z]+)(\d+)", atom2)
    atom2_label = match.group(1)
    atom2_id = match.group(2)

    if atom1_label < atom2_label:
        bond_key = atom1_label + "-" + atom2_label
    else:
        bond_key = atom2_label + "-" + atom1_label

    bond_length_dict = avg_bond_lengths[bond_key]
    closest = 'single'
    for (order, length) in bond_length_dict.items():
        length_diff = abs(bond_length - float(length))
        if length_diff < abs(bond_length_dict[closest] - bond_length):
            closest = order

    return (atom1_label, atom1_id, atom2_label, atom2_id, closest)

# table of type: 
# atom1 label, atom1 id, atom2 label, atom2 id, bond order
def table_to_gml(table):
    gml_bond_char = {'single': '-', 'double':'=', 'triple':'#'}
    lines = ["graph ["]
    nodes = []
    for (atom1_label, atom1_id, atom2_label, atom2_id, order) in table:
        if atom1_id not in nodes:
            lines.append(f"\tnode [ id {atom1_id} label \"{atom1_label}\" ]")
            nodes.append(atom1_id)
        if atom2_id not in nodes:
            lines.append(f"\tnode [ id {atom2_id} label \"{atom2_label}\" ]")
            nodes.append(atom2_id)

        lines.append(f"\tedge [ source {atom1_id} target {atom2_id} label \"{gml_bond_char[order]}\" ]")

    lines.append("]")
    return "\n".join(lines)

def write_gml_file(bond_table, filename="unnamed"):
    classification = [classify_bond(x,y) for (x,y) in zip(bond_table['A-B'], bond_table['distance_calc'])]
    with open(f"{filename}.gml", 'w') as file:
        file.write(table_to_gml(classification))

def read_allfrags(args, allfrags_dir="."):
    Path(f"./iso_fragments").mkdir(exist_ok=True)
    for f in Path(f"./iso_fragments").glob("*"):
        if f.is_file():
            f.unlink()

    Path(f"./pair_fragments").mkdir(exist_ok=True)
    for f in Path(f"./pair_fragments").glob("*"):
        if f.is_file():
            f.unlink()

    with open(f"{allfrags_dir}/allfragments") as f:
        line = f.readline() # discard header
        while line:
            line = f.readline()
            tokens = line.split()
            if len(tokens) != 6:
                continue
            [dire, frag_type, _, _,  _, _] = tokens
            if frag_type=='isomer':
                pt = PrintTab(args, f"{allfrags_dir}/{dire}/isomer.xyz")
                write_gml_file(pt.bond_table, f"./iso_fragments/{dire}")
                line = f.readline() # skip next line
            elif frag_type=='fragmentpair':
                [frag1_dir, _, _, _] = f.readline().split()
                [frag2_dir, _, _, _] = f.readline().split()

                try:
                    pt1 = PrintTab(args, f"{allfrags_dir}/{frag1_dir}/fragment.xyz")
                    write_gml_file(pt1.bond_table, f"./pair_fragments/{frag1_dir}")
                except:
                    pass
                
                try:
                    pt2 = PrintTab(args, f"{allfrags_dir}/{frag2_dir}/fragment.xyz")
                    write_gml_file(pt2.bond_table, f"./pair_fragments/{frag2_dir}")
                except:
                    pass

            

def main():
    args = parse_args()
    
    #for windows console
    sys.stdout.reconfigure(encoding='utf-8')  

    read_allfrags(args, allfrags_dir=args.allfrags_dir)
    

if __name__=="__main__":
    main()