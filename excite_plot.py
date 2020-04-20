import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.pyplot as plt

font = {'size':22}
mpl.rc('font',**font)

infile = sys.argv[1]


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

df = loadExciteDf(infile)
df = df.assign(E_nm=hc/df['E_eV'])
df = df.sort_values(by=['Osc_Str'])
print(df)

x_arr = np.linspace(350,800,n_grid)
spectrum = np.zeros(n_grid)
Fe_spectrum = np.zeros(n_grid)

for index, row in df.iterrows():
    spectrum += gaussian(x_arr, row['E_nm'], row['Osc_Str']**2, sig2)


plt.plot(x_arr, spectrum, label='PhencirFe-O-PhencirFe')
plt.ylabel('Absorbance')
plt.xlabel('Wavelength (nm)')
plt.legend()
plt.show()
#plt.savefig("nan_spectrum.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
