import Qdata
import ChemData

import pandas as pd
import sys
import os
from collections import Counter

"""
Calculate dG_formation of molecule using a solvated calculation at optimum geometry and a frequency calculation
Assumes Q-Chem output file format

When calling this script directly:
Input 1: directory containing all the necessary Q-Chem jobs
Input 2: name of output csv file
"""

class Gform(object):
  def __init__(self, in_directory, outfile_name, therm_datafile):
    print(therm_datafile)
    therm_col = ['mol_tag','E','dHf0','ddH-ddH_ss','H_zpe','S-S_ss','dGfaq']
    self.Gf_df = pd.DataFrame(columns=therm_col)
    self.therm_data = pd.read_csv(therm_datafile, index_col=0, dtype={'S':float}) #Contains thermo ref info 
    print(self.therm_data)
    self.indir = in_directory
    self.outfile_name = outfile_name
    #Conversion factors and constants
    self.T = 298.15 #The assumed temperature
    self.Ht2kC = 627.509469 #The conversion factor for Hartree to kCal/mol
    self.vib_scale = 0.9512 #An empirical parameter for fitting DFT freq to experiment
    self.kCal2kJ = 4.184 
    self.Htrpv = 2.372
    #SSmPa is a dictionary with molecules per atom in standard state
    self.SSmPa = {'H':0.5,'C':1.0,'N':0.5,'O':0.5,'F':0.5,'Ni':1.0,'S':1.0,'Fe':1.0}

#Parses qchem output files and loads a Qdata instance with appropriate data
  def dG(self, fq_file=None, pcm_file=None):
    print("Performing PCET Calc: ", fq_file, pcm_file)
    qdat = Qdata.Qdata()
    qdat.readFile(fq_file)
    #qdat.H = qdat_fq.H
    #qdat.Hzpe = qdat_fq.Hzpe
    #qdat.Hvib = qdat_fq.Hvib
    #qdat.S = qdat_fq.S
    if qdat.H == None: raise ValueError("Couldn't find necessary thermo info")
    if pcm_file != None: 
        qdat_solv = Qdata.Qdata()
        qdat_solv.readFile(pcm_file)
        qdat.E = qdat_solv.E
    return self.calcG(qdat) #overwrites qdat.E with E including solv effects

  def calcG(self, qdat):
    Eatoms, dHfatoms, ddH_ss, S_ss = self.thermSS(qdat)
    ddH = qdat.Hvib*self.vib_scale - qdat.Hzpe*self.vib_scale + self.Htrpv
    dHAtomization = -1.0*qdat.E*self.Ht2kC + Eatoms*self.Ht2kC - qdat.Hzpe*self.vib_scale #solvation energy includes total E
#    dHAtomization = -1.0*qdat.Esolv*self.Ht2kC + Eatoms*self.Ht2kC - qdat.Hzpe*self.vib_scale
    dHf0 = dHfatoms - dHAtomization
    dHf298 = dHf0 + ddH - ddH_ss
    dGfg = dHf298 - self.T*(qdat.S - S_ss)
    print('Eatoms: ',Eatoms)
    print('dHfatoms: ', dHfatoms)
    print('dHf0: ', dHf0)
    print('ddH_ss: ', ddH_ss)
    print('S_ss: ', S_ss)
    print('E: ', qdat.E)
    print('H_zpe: ', qdat.Hzpe)
    print('ddH: ',ddH)
    print('S: ',qdat.S)
    print('dHatomization: ',dHAtomization)
    dGfaq = dGfg + 1.89 #+ qdat.Gpcm #qdat.loadPCM overwrites qdat.E w/ E including solv effects
    return [qdat.E,dHf0,ddH-ddH_ss,qdat.Hzpe,qdat.S-S_ss,dGfaq]

#Main workhorse; iterate through target directory and scrape out thermodynamic information
  def batchDG(self):
    for filename in os.listdir(self.indir):
      if filename.split('.')[-1] == 'out' and 'fq' in filename: #ignore all non-fq output files
        job = filename.split('_')[1]
  
        #if job == 'sfq': #sfq is a solvated frequency calculation; this has all the info we need
        #  therm_list = self.dG(self.indir+"/"+filename,pcm_file=(self.indir+"/"+filename))
        #elif job in ['opfqSo', 'sopfqSo', 'optspfq']: #opt-freq-solventSP 
        #  therm_list = self.dG(self.indir+'/'+filename)
  
        #freq run separately from solvation; look for corresponding PCM file
        chem_id = self.chemID(filename)
        for pcm_file in os.listdir(self.indir):
          if ('sp' or 'Sm8' or 'Sm12' or 'sop') in pcm_file :
            pcm_id = self.chemID(pcm_file)
            if chem_id == pcm_id:
              therm_list = self.dG(self.indir+"/"+filename, pcm_file=self.indir+"/"+pcm_file)
              break

      else: 
        print("Just not sure where this job is going...")
        continue #whatever the file is, it doesn't follow a recognizeable naming structure

      print(filename, therm_list[-1]) #last item in therm_list is dGfaq, the important quantity
      self.Gf_df.loc[len(self.Gf_df)] = [filename]+therm_list
    print()
    print(self.Gf_df)
    self.Gf_df.to_csv(self.outfile_name)

#the first part is whatever I named the molecule, the last is its spin-charge tag
#this function just takes those 2 pieces because they should be enough to uniquely indentify the system
  def chemID(self, qoutname):
    splile = qoutname.split('.')[0]
    splile = splile.split('_')
    return splile[0] + splile[-1]

  def thermSS(self, Qdata):
    """
    calculates thermodynamic data for the steady state
    output units: 
    Eatoms: Hartree 
    dHfatoms, ddH_ss, S_ss: kCal/mol
    """
    dHfatoms = 0.0
    S_ss = 0.0
    ddH_ss = 0.0
    Eatoms = 0.0
    atom_count = Counter(Qdata.atoms)
    print(atom_count)
    for key_atom in atom_count:
      dHfatoms += atom_count[key_atom] * self.therm_data['H0'][key_atom]/self.kCal2kJ
      ddH_ss += atom_count[key_atom] * self.therm_data['dHss'][key_atom]*\
                self.SSmPa[key_atom]/self.kCal2kJ
      print(type(atom_count[key_atom]),type(self.therm_data['S'][key_atom]),type(self.SSmPa[key_atom]))
      print(self.therm_data['S'][key_atom])
      print(type(self.therm_data['S'][key_atom]))
      print(key_atom)
      print(self.therm_data['S']['O'])
      S_ss += atom_count[key_atom] * self.therm_data['S'][key_atom]*self.SSmPa[key_atom]/self.kCal2kJ
      Eatoms += atom_count[key_atom] * self.therm_data['E'][key_atom]
    return Eatoms, dHfatoms, ddH_ss ,S_ss

if __name__ == "__main__":
  in_dir = sys.argv[1] #directory where all of the jobs are stored
  out_df_name = sys.argv[2] #outputs pandas dataframe to csv
  #tdat = os.path.join(sys.path[0],"tdat_tpssh.csv")
  #tdat = os.path.join(sys.path[0],"tdat_pbe.csv")
  tdat = os.path.join(sys.path[0],"tdat_b3lyp.csv")
  gform = Gform(in_dir,out_df_name,tdat)
  gform.batchDG() 




