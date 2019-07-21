import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.pyplot as plt

font = {'size':22}
mpl.rc('font',**font)

#infile = sys.argv[1]

def loadExciteDf(infile):
    strength_list = []
    energy_list = []
    with open(infile, 'r') as f:
        for line in f.readlines():
            if 'excitation energy' in line:
                energy_list.append(float(line.split()[-1]))
            if 'Strength' in line:
                strength_list.append(float(line.split()[-1]))

    excit_dict = {'Osc_Str':strength_list, 'E_eV':energy_list}
    df = pd.DataFrame(excit_dict)
    return df

def gaussian(x, mu, norm, sig2):
    return norm/((2.*np.pi*sig2)**0.5)*np.exp(-1.*(x-mu)**2/(2*sig2))


hc = 1240.
n_grid = 2500
sig2 = 20.

nan_df = loadExciteDf('nan_tddft_a0m1_states.txt')
nanFeCl_df = loadExciteDf('nanFeCl_tddft_a0m4_states.txt')

nan_df = nan_df.assign(E_nm=hc/nan_df['E_eV'])
nanFeCl_df = nanFeCl_df.assign(E_nm=hc/nanFeCl_df['E_eV'])

nan_df = nan_df.sort_values(by=['Osc_Str'])
nanFeCl_df = nanFeCl_df.sort_values(by=['Osc_Str'])

print(nan_df)
print(nanFeCl_df)

x_arr = np.linspace(150,500,n_grid)
spectrum = np.zeros(n_grid)
Fe_spectrum = np.zeros(n_grid)

for index, row in nan_df.iterrows():
    spectrum += gaussian(x_arr, row['E_nm'], row['Osc_Str']**2, sig2)

for index, row in nanFeCl_df.iterrows():
    Fe_spectrum += gaussian(x_arr, row['E_nm'], 22.*row['Osc_Str']**2, sig2)

plt.plot(x_arr, spectrum, label='Phencir')
plt.plot(x_arr, Fe_spectrum, label='PhencirFe')
plt.ylabel('Absorbance')
plt.xlabel('Wavelength (nm)')
plt.legend()
plt.show()
#plt.savefig("nan_spectrum.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
