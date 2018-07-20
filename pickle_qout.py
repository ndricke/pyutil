# input: file or directory
# pickles Q-Chem .out files

import sys
import os
import Qdata
import pickle
import dill

def read_data(qchem_file):
    if qchem_file.split('.')[-1] == 'out':
        qdata = Qdata.Qdata()
        qdata.readFile(qchem_file)
        return qdata


in_filedir = sys.argv[1]
pickle_fname = sys.argv[2]

save_data = []
for subdir, dirs, files in os.walk(in_filedir):
    for qfile in files:
        print(os.path.join(subdir,qfile))
        if qfile.split('.')[-1] == 'out':
            qdata = read_data(os.path.join(subdir,qfile))
            save_data.append(qdata)
        
dill.dump(save_data, open(pickle_fname, 'wb'))
