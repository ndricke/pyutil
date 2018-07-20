import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#class QchemData

def sectionGrab(key_start,key_end,in_file,start_cut=0,end_cut=0):
 y_read = 0
 grab_list = []
 f = open(in_file,'r')
 for line in f:
  if key_start in line:
   y_read = 1
  elif key_end in line:
   y_read = 0
  elif y_read == 1:
   grab_list.append(line.strip())
 f.close()
 if end_cut == 0:
  return grab_list[start_cut:]
 else:
  return grab_list[start_cut:end_cut]

def chargeGrab(in_file):
 grab_start = "Ground-State Mulliken Net Atomic Charges"
 grab_end = "Sum of atomic charges"
 charge_raw = sectionGrab(grab_start,grab_end,in_file,3,-1)
 charge = np.zeros(len(charge_raw))
 atom = []
 i = 0
 for line in charge_raw:
  spline = line.split()
  charge[i] = float(spline[2])
  atom.append(spline[1])
  i += 1
 return charge, atom

def beckeGrab(in_file):
 grab_start = "CDFT Becke Populations"
 grab_end = "32       H"
 charge_raw = sectionGrab(grab_start,grab_end,in_file,2,-1)
 print(charge_raw)
 charge = np.zeros(len(charge_raw))
 atom = []
 i = 0
 for line in charge_raw:
  spline = line.split()
  charge[i] = float(spline[2])
  atom.append(spline[1])
  i += 1
 return charge, atom

def coordGrab(in_file):
 grab_start = 'Standard Nuclear Orientation (Angstroms)'
 grab_end = 'Nuclear Repulsion Energy'
 raw_coord = sectionGrab(grab_start,grab_end,in_file,2,-1)
 coord = np.zeros((len(raw_coord),3))
 i = 0
 for line in raw_coord:
  spline = line.split()
  coord[i,:] = [float(j) for j in spline[2:]]
  i += 1
 return coord

atom_dict = {'C':'k','N':'b','F':'r','H':'w'}
size_scale = 2000

qfile = sys.argv[1]

#a_charge, a_atom = chargeGrab(qfile)
a_charge, a_atom = beckeGrab(qfile)
a_coord = coordGrab(qfile)

print(a_charge)
print(a_atom)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

a_atom_color = []
for atom in a_atom:
 a_atom_color.append(atom_dict[atom])

a_mark = []
for charge in a_charge:
 if np.sign(charge) >= 0:
  a_mark.append('o')
 else:
  a_mark.append('v')

for i in range(len(a_charge)):
 ax.scatter(a_coord[i,0], a_coord[i,1],a_coord[i,2],s=np.abs(a_charge[i])*size_scale, \
            c=a_atom_color[i],marker=a_mark[i] )

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
#ax.title(mol_name)

#plt.savefig(out_dir + '/' + mol_name + '_' + diff_name + '.png')
plt.show()
plt.clf()

























