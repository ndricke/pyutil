import os

dir_choice = "inp"
dir_new = "out"

for filename in os.listdir(dir_choice):
 if filename[-4:] == ".txt":
  y_read = 0
  fpath = '/'.join([dir_choice, filename])
  outname = '.'.join([filename[:-4],'xyz'])
  wpath = '/'.join([dir_new, outname])
  print(fpath)
  k = sum(1 for line in open(fpath,'r'))
  f = open(fpath,'r')
  newfile = open(wpath,'w')
  newfile.write(str(k))
  newfile.write('\n')
  newfile.write('\n')
  for line in f:
   ban = line.split()
   nline = ' '.join(ban[1:])
   newfile.write(nline + '\n')
  f.close()
  newfile.close()
   
