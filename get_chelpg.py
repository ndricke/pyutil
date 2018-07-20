# input: pickled directory or file
# creates DataFrame containing file names and identities/CHELPG charges
    # of any metal atoms in catalyst

import numpy as np
import pandas as pd
import sys
import os
import Qdata
import pickle

in_pickle = sys.argv[1]
save_data = pickle.load(open(in_pickle, "rb"))
metals = ["Fe", "Mn", "Co", "Cr"] # for modularity, make this a parameter instead
filenames, atoms, chelpgs = [], [], []

print(in_pickle)

for qdata in save_data:
    try:
        chelpg = [float(item) for item in qdata.chelpg]
        chelpg_df = pd.DataFrame({'Atom':qdata.atoms,'CHELPG':chelpg})
        for index, row in chelpg_df.iterrows():
            for metal in metals:
                if row['Atom'] == metal:
                    filenames.append(qdata.filename)
                    atoms.append(row['Atom'])
                    chelpgs.append(row['CHELPG'])
    except:
        continue
chelpg_metals = pd.DataFrame.from_items([('File_Name', filenames), ('Atom', atoms), ('CHELPG', chelpgs)])
print(chelpg_metals)
