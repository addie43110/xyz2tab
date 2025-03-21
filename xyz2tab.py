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
    bond_table = pt.bond_table
    if bond_table:
        classification = [classify_bond(x,y) for (x,y) in zip(bond_table['A-B'], bond_table['distance_calc'])]
        gml_string = table_to_gml(classification)
    else: #there is only a single atom, so could not make any bond information
        element = pt.xyz_df[['element']][0]
        gml_string = f"graph [\n\tnode [ id 0 label {element} ]\n]"
    try:
        mod_graph = mod.Graph.fromGMLString(gml_string)
    except mod.libpymod.InputError:
        print(f"Error trying to write {filename}.gml. Likely graph is not connected. Skipping.")
        return None
    g = Graph(modGraph=mod_graph)
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
    match = re.match(r"^(\w*)([pf]\d+)(?:p\d+)(?:f\d+)?$", directory)
    look_in_iso = True if match.group(2)[0]=='p' else False
    fragment_dir = "iso_fragments" if look_in_iso else "pair_fragments"
    parent_filename = match.group(1)+match.group(2)
    try:
        parent_graph = mod.Graph.fromGMLFile(f"./{fragment_dir}/{parent_filename}.gml")
    except:
        print("Parent fragment not found. Omitting rules with direct children.")
        return (None, None)
    parent_graph = Graph(modGraph=parent_graph)
    parent_name = parent_filename
    return (parent_graph, parent_name)


def read_allfrags(args, allfrags_dir=".", initial_pname="unnamed"):
    make_exist_dir("./iso_fragments")
    make_exist_dir("./pair_fragments")
    make_exist_dir("./rules")

    parent_name=initial_pname
    pt_start = PrintTab(args, f"{allfrags_dir}/in.xyz")
    parent_graph = write_gml_file(pt_start, f"./iso_fragments/{initial_pname}")
    update_parent = False

    with open(f"{allfrags_dir}/allfragments") as f:
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
                path_to_fragment = f"{allfrags_dir}/{dire}/isomer.xyz"
                where_to_write_gml = f"./iso_fragments/{dire}"
                where_to_write_rule = f"./rules/{parent_name}_{dire}.gml"
                read_fragment(args, path_to_fragment, where_to_write_gml, where_to_write_rule, parent_graph, parent_name)
                f.readline() # skip next line

            elif frag_type=='fragmentpair':
                [frag1_dir, _, _, _] = f.readline().split()
                [frag2_dir, _, _, _] = f.readline().split()

                ptf = [f"{allfrags_dir}/{d}/fragment.xyz" for d in [frag1_dir, frag2_dir]]
                wtwg = [f"./pair_fragments/{d}" for d in [frag1_dir, frag2_dir]]
                wtwr = [f"./rules/{parent_name}_{d}.gml" for d in [frag1_dir, frag2_dir]]

                read_fragment(args, ptf[0], wtwg[0], wtwr[0], parent_graph, parent_name)
                read_fragment(args, ptf[1], wtwg[1], wtwr[1], parent_graph, parent_name)

            elif frag_type=='fragment_type':
                update_parent = True

def read_fragment(args, path_to_fragment, where_to_write_gml, where_to_write_rule, parent_graph, parent_name):
    """ try:
        pt = PrintTab(args, path_to_fragment)
    except Exception as e:
        print(f"Encountered malformed .xyz file at {path_to_fragment}. Skipping.")
        print(e)
        return """
    pt = PrintTab(args, path_to_fragment)
    child_graph = write_gml_file(pt, where_to_write_gml)
    if child_graph and parent_graph:
        rule_gml_string = Reaction(leftGraph=parent_graph, rightGraph=child_graph, name=where_to_write_rule[8:-4]).to_ruleGML_string()
        write_gml_string(rule_gml_string, where_to_write_rule)

def main():
    args = parse_args()
    
    #for windows console
    sys.stdout.reconfigure(encoding='utf-8')  

    read_allfrags(args, allfrags_dir=args.allfrags_dir[:-1], initial_pname=args.name)
    

if __name__=="__main__":
    main()