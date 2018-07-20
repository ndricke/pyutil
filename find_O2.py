# input: pickled directory or file
# identifies exactly two oxygen atoms as those in dioxygen; warning otherwise
# determines catalyst atom to which O2 is bound
# creates DataFrame containing file names and indices for these three
    # atoms of interest in each catalyst, or "NaN" if not applicable
# DataFrame indices refer to atom_coords list from pickle, NOT list
    # of atoms in Q-Chem .out file. Q-Chem index = atom_coords index + 1

import numpy as np
import pandas as pd
import sys
import os
import Qdata
import dill

in_pickle = sys.argv[1]
save_data = dill.load(open(in_pickle, "rb"))

filenames, oxy1_list, oxy2_list, cat_list, energy_list = [], [], [], [], []
M_O2_bond_lengths = []
# used to create DataFrame

def original_index(O_index):
    # converts O_coords index to atom_coords index; precondition: O_index >= 0
    count = 0
    for i in range(len(atom_coords)):
        if atom_coords[i].split(' ')[0] == 'O':
            count += 1
        if count == O_index + 1:
            return i
    print(O_index + " is not a valid O_coords index.")

def get_coords(input_string):
    # input_string of the form "<atom> <x-coord> <y-coord> <z-coord>"
    return np.array([float(input_string.split(' ')[1]), float(input_string.split(' ')[2]), float(input_string.split(' ')[3])])

def distance(string1, string2):
    difference = get_coords(string1) - get_coords(string2)
    return abs(np.linalg.norm(difference))

for qdata in save_data:
    filenames.append(qdata.filename)
    try:
        energy_list.append(qdata.E)
    except:
        energy_list.append(None)
    O2_found, multiple_O2 = False, False

    try:
        atom_coords = qdata.listCoord()
    except:
        print("No coordinates in: ", qdata.filename)
        oxy1_list.append(None)
        oxy2_list.append(None)
        M_O2_bond_lengths.append(None)
        cat_list.append(None)
        continue
    O_coords = []
    for coord in atom_coords:
        if coord.split(' ')[0] == 'O':
            O_coords.append(coord)
    if len(O_coords) < 2:
        #print("Fewer than two oxygens in " + qdata.filename)
        oxy1_list.append(None)
        oxy2_list.append(None)
        M_O2_bond_lengths.append(None)
        cat_list.append(None)
        continue
    for index1 in range(len(O_coords) - 1):
        for index2 in range(index1 + 1, len(O_coords)):
            if distance(O_coords[index1], O_coords[index2]) < 1.6:
                if O2_found == True:
                    multiple_O2 = True
                    break
                    # doesn't exit index1 for loop but still reduces unnecessary operations
                else:
                    O2_found = True
                    oxy1, oxy2 = index1, index2 # O_coords index values for oxygens identified as O2
    if not O2_found:
        #print("No O2 detected in " + qdata.filename)
        oxy1_list.append(None)
        oxy2_list.append(None)
        M_O2_bond_lengths.append(None)
        cat_list.append(None)
        continue
    if multiple_O2:
        #print("Multiple O2 detected in " + qdata.filename)
        oxy1_list.append(None)
        oxy2_list.append(None)
        M_O2_bond_lengths.append(None)
        cat_list.append(None)
        continue

    shortest_dist, catalyst_index, oxygen_index = sys.float_info.max, sys.float_info.max, sys.float_info.max
    # catalyst_index and oxygen_index refer to atom_coords, not O_coords
    oxy1_orig, oxy2_orig = original_index(oxy1), original_index(oxy2)
    oxy1_list.append(oxy1_orig)
    oxy2_list.append(oxy2_orig)
    # ensures oxy1_list, oxy2_list have correct length; last entry may be overwritten
    for atom_index in range(len(atom_coords)):
        if atom_index == oxy1_orig or atom_index == oxy2_orig:
            continue
        if distance(atom_coords[oxy1_orig], atom_coords[atom_index]) < shortest_dist:
            shortest_dist = distance(atom_coords[oxy1_orig], atom_coords[atom_index])
            catalyst_index = atom_index
            oxygen_index = oxy1_orig
            oxy1_list[-1] = oxy1_orig
            oxy2_list[-1] = oxy2_orig
        if distance(atom_coords[oxy2_orig], atom_coords[atom_index]) < shortest_dist:
            shortest_dist = distance(atom_coords[oxy2_orig], atom_coords[atom_index])
            catalyst_index = atom_index
            oxygen_index = oxy2_orig
            oxy1_list[-1] = oxy2_orig
            oxy2_list[-1] = oxy1_orig
            # oxy1_list contains the O bound to the catalyst, oxy2_list the other O
    M_O2_bond_lengths.append(shortest_dist)
    cat_list.append(catalyst_index)

O2_df = pd.DataFrame.from_items([('File_Name', filenames), ('Energy', energy_list), ('Active_Site', cat_list), ('Oxygen_1', oxy1_list), 
                                 ('Oxygen_2', oxy2_list), ('M-O2_Bond_Length', M_O2_bond_lengths)])

print(O2_df)
O2_df.to_csv('pyrpyc.csv')
O2_df.to_pickle('pyrpyc_df.p')



























