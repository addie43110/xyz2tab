#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#xyz2tab

import sys
import re
import mod
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

from graph import Graph
from reaction import Reaction
from lookup_tables import avg_bond_lengths
from helpers import parse_args
from print_tables import PrintTab
from prettify import red, warn, green, blue
# from create_plots import CreatePlot

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

# returns Graph object of the gml file written
def write_gml_file(pt, filename="unnamed") -> Graph:
    if pt.has_bond_table:
        bond_table = pt.bond_table
        classification = [classify_bond(x,y) for (x,y) in zip(bond_table['A-B'], bond_table['distance_calc'])]
        gml_string = table_to_gml(classification)
    else: #there is only a single atom, so could not make any bond information
        element = pt.xyz_df.iloc[0]['element']
        gml_string = f"graph [\n\tnode [ id 0 label {element} ]\n]"

    g = Graph(gml_string)

    try:
        mod_graph = mod.Graph.fromGMLString(gml_string)
        g = Graph(mod_graph)
    except mod.libpymod.InputError:
        print(f"Error trying to write {filename}.gml. Likely graph is not connected or no edges found.")

    with open(f"{filename}.gml", 'w') as file:
        file.write(gml_string)
    return g

def write_gml_string(gml_string, filename="./unamed.gml"):
    with open(f"{filename}", 'w') as file:
        file.write(gml_string)

def make_exist_dir(dir_path):
    Path(dir_path).mkdir(exist_ok=True)
    for f in Path(dir_path).glob("*"):
        if f.is_file():
            f.unlink()

def updateParent(directory):
    try:
        parent_graph = mod.Graph.fromGMLFile(f"./all_fragments/{parent_filename}.gml")
    except:
        print("Parent fragment not found. Omitting rules with direct children.")
        return (None, None)
    parent_graph = Graph(modGraph=parent_graph)
    parent_name = parent_filename
    return (parent_graph, parent_name)

def read_peakfrags(qcxsm2_dir):
    peak_dict = {}
    with open(f"{qcxsm2_dir}/allpeaks.dat") as f:
        f.readline()
        line = f.readline() # discard two headers
        while line:
            line = f.readline()
            tokens = line.split()
            if len(tokens) != 3:
                continue
            (frag_name, _, intensity) = tokens
            peak_dict[frag_name] = float(intensity)
    return peak_dict

def read_allfrags(args, qcxsm2_dir=".", initial_pname="unnamed"):
    peak_dict = read_peakfrags(qcxsm2_dir)

    make_exist_dir("./all_fragments")
    make_exist_dir("./peak_fragments")
    make_exist_dir("./rules")

    parent_name=initial_pname
    pt_start = PrintTab(args, f"{qcxsm2_dir}/in.xyz")
    parent_graph = write_gml_file(pt_start, f"./all_fragments/{initial_pname}")
    update_parent = False

    with open(f"{qcxsm2_dir}/allfragments") as f:
        line = f.readline() # discard header
        while line:
            line = f.readline()
            tokens = line.split()
            if len(tokens) != 6:
                continue
            [dire, frag_type, _, _,  _, _] = tokens
            if update_parent:
                    (parent_graph, parent_name) = updateParent(dire)
                    update_parent = False
            if frag_type=='isomer':
                path_to_fragment = f"{qcxsm2_dir}/{dire}/isomer.xyz"
                read_fragment(args, path_to_fragment, dire, parent_graph, parent_name, peak_dict)
                f.readline() # skip next line

            elif frag_type=='fragmentpair':
                [frag1_dir, _, _, _] = f.readline().split()
                [frag2_dir, _, _, _] = f.readline().split()

                ptf = [f"{qcxsm2_dir}/{d}/fragment.xyz" for d in [frag1_dir, frag2_dir]]
                read_fragment(args, ptf[0], frag1_dir, parent_graph, parent_name, peak_dict)
                read_fragment(args, ptf[1], frag2_dir, parent_graph, parent_name, peak_dict)

            elif frag_type=='fragment_type':
                update_parent = True

def read_fragment(args, path_to_fragment, frag_name, parent_graph, parent_name, peak_dict):
    pt = PrintTab(args, path_to_fragment)
    child_graph = write_gml_file(pt, f"./all_fragments/{frag_name}")
    if peak_dict[frag_name] >= 1:
        write_gml_file(pt, f"./peak_fragments/{frag_name}")
    if child_graph and parent_graph and peak_dict[frag_name] >= 1:
        rule_gml_string = Reaction(leftGraph=parent_graph, rightGraph=child_graph, name=f"{parent_name}_{frag_name}").to_ruleGML_string()
        write_gml_string(rule_gml_string, f"./rules/{parent_name}_{frag_name}")

def main():
    args = parse_args()
    
    #for windows console
    sys.stdout.reconfigure(encoding='utf-8')  

    read_allfrags(args, qcxsm2_dir=args.allfrags_dir[:-1], initial_pname=args.name)
    

if __name__=="__main__":
    main()