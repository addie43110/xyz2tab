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
from openbabel import openbabel

from graph import Graph
from reaction import Reaction
from lookup_tables import avg_bond_lengths
from helpers import parse_args
from print_tables import PrintTab
from prettify import red, warn, green, blue
# from create_plots import CreatePlot

#pd.set_option("display.max_rows", None, "display.max_columns", None)

# Classifies bond order based on bond length
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
    closest = list(bond_length_dict.keys())[0]
    for (order, length) in bond_length_dict.items():
        length_diff = abs(bond_length - float(length))
        if length_diff < abs(bond_length_dict[closest] - bond_length):
            closest = order

    return (atom1_label, atom1_id, atom2_label, atom2_id, closest)

# table of type: 
# atom1 label, atom1 id, atom2 label, atom2 id, bond order
def table_to_gml(table, pt=None):
    gml_bond_char = {'single': '-', 'double':'=', 'triple':'#', 'aromatic':':'}
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

    if pt.num_atoms > len(nodes):
        lines.append("")
        for atom_idx in pt.xyz_df['atom1_idx']:
            m = re.match(r"([a-zA-Z]+)(\d+)", atom_idx)
            atom_label = m.group(1)
            atom_id = m.group(2)
            if atom_id not in nodes:
                lines.append(f"\tnode [ id {atom_id} label \"{atom_label}\" ]")
    lines.append("]")
    return "\n".join(lines)

# general write gml string (graph or rule)
def write_gml_string(gml_string, filename="./unamed.gml"):
    with open(f"{filename}", 'w') as file:
        file.write(gml_string)

def make_exist_dir(dir_path):
    Path(dir_path).mkdir(exist_ok=True)
    for f in Path(dir_path).glob("*"):
        if f.is_file():
            f.unlink()

""" directory is a string given by QCxMS2 to identify fragments.
    For example: p1p15 is the 15th product of fragmentation of p1
                 p11f1 is the first fragment of p11 (p11 is a pair of fragments)
    Regex is used to determine the parent product.
    In the first example, p1 is the parent product.
    In the second example, the original molecule is the parent. """
def updateParent(directory):
    match = re.match(r"^(\w*)([pf]\d+)(?:p\d+)(?:f\d+)?$", directory)
    parent_filename = match.group(1)+match.group(2)
    try:
        parent_graph = mod.Graph.fromGMLFile(f"./all_fragments/{parent_filename}.gml")
    except:
        print(f"Parent fragment {parent_filename} not found. Omitting rules with direct children.")
        # print(f"Printing contents of {parent_dir}...")
        # print(s)
        return (None, None)
    parent_graph = Graph(parent_graph)
    return (parent_graph, parent_filename)

""" Opens the file 'allpeaks.dat' found in the root directory of the 
    QCxMS2 data. Creates a dictionary where keys are the fragment name
    and values are the intensity of the peak. """
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
    parent_gml = pt_to_gml(args, f"{qcxsm2_dir}/in.xyz")
    xyz_to_gml(f"{qcxsm2_dir}/in.xyz")
    #parent_m = Chem.MolFromXYZFile(f"{qcxsm2_dir}/in.xyz")
    #print(f"rdkit mol type: {type(parent_m)}")
    write_gml_string(parent_gml, f"./all_fragments/{parent_name}.gml")
    parent_graph = Graph(parent_gml)
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
                path_to_pair = f"{qcxsm2_dir}/{dire}/pair.xyz"
                read_fragment(args, path_to_pair, dire, parent_graph, parent_name, peak_dict)
                f.readline()
                f.readline() # skip next two lines

            elif frag_type=='fragment_type':
                update_parent = True

def xyz_to_gml(path_to_xyz):
    match = re.match(r"^(\S+/)*(\w+.xyz)$", path_to_xyz)
    parent_dirs = match.group(1)
    filename = match.group(2)[:-4]

    conv_obj = openbabel.OBConversion()
    conv_obj.SetInFormat("xyz")
    mol = openbabel.OBMol()
    conv_obj.ReadFile(mol, path_to_xyz)

    #for bond in openbabel.OBMolBondIter(mol):
    #    print(f"start: {bond.GetBeginAtomIdx()}, end: {bond.GetEndAtomIdx()}, length: {bond.GetLength()}, order: {bond.GetBondOrder()}")

    charge_model = openbabel.OBChargeModel.FindType("eem")
    print(f"\ncharge computed?: {charge_model.ComputeCharges(mol)}")
    print(f"formal charges: {charge_model.GetFormalCharges()}")
    

def pt_to_gml(args, path_to_fragment):
    pt = PrintTab(args, path_to_fragment)
    if pt.has_bond_table:
        bond_table = pt.bond_table
        classification = [classify_bond(x,y) for (x,y) in zip(bond_table['A-B'], bond_table['distance_calc'])]
        gml_string = table_to_gml(classification, pt)
    else: #there is only one or two atoms which have no bonds
        lines = []
        lines.append("graph [")
        for atom_idx in pt.xyz_df['atom1_idx']:
            m = re.match(r"([a-zA-Z]+)(\d+)", atom_idx)
            atom_label = m.group(1)
            atom_id = m.group(2)
            lines.append(f"\tnode [ id {atom_id} label \"{atom_label}\" ]")
        lines.append("]")
        gml_string = "\n".join(lines)
    return gml_string
    

def read_fragment(args, path_to_fragment, frag_name, parent_graph, parent_name, peak_dict):
    gml_string = pt_to_gml(args, path_to_fragment)
    xyz_to_gml(path_to_fragment)
    g = Graph(gml_string)
    ccps = g.connected_components

    is_peak = False
    if len(ccps)==1: #isomer
        write_gml_string(gml_string, f"./all_fragments/{frag_name}.gml")
        if peak_dict[frag_name] >= 1:
            write_gml_string(gml_string, f"./peak_fragments/{frag_name}.gml")
            is_peak = True
    elif len(ccps)==2: # pair fragments
        frag_suffixes = ["f1", "f2"]
        for fs, cmp in zip(frag_suffixes, ccps):
            write_gml_string(str(cmp), f"./all_fragments/{frag_name}{fs}.gml")
            if peak_dict[f"{frag_name}{fs}"] >= 1:
                write_gml_string(str(cmp), f"./peak_fragments/{frag_name}{fs}.gml")
                is_peak = True
    else:
        print(red(f"ERROR: file {frag_name} contains more than 2 fragments. Not writing gml files."))

    # test: maybe need to write all rules regardless of is_peak or parent.
    #if parent_graph and is_peak:
    rule_gml_string = Reaction(educts=[parent_graph], products=ccps, name=f"{parent_name}!!{frag_name}").to_ruleGML_string()
    write_gml_string(rule_gml_string, f"./rules/{parent_name}_{frag_name}.gml")

def main():
    args = parse_args()
    
    #for windows console
    sys.stdout.reconfigure(encoding='utf-8')  

    read_allfrags(args, qcxsm2_dir=args.allfrags_dir[:-1], initial_pname=args.name)
    

if __name__=="__main__":
    main()