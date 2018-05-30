import sys
import numpy as np
import matplotlib.pyplot as plt
import os

import Qdata

def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def orbitals2Density(mo_energies, dom, sig):
    density = np.zeros(dom_num)
    for mu in mo_energies:
        density += gaussian(dom, mu, sig)
    return density

def genDensity(infile, domain, sigma=0.005):
    qdata = Qdata.Qdata()
    qdata.readFile(infile)

    alpha_occ = [float(item) for item in qdata.alpha_occ[1:]]
    beta_occ = [float(item) for item in qdata.beta_occ[1:]]
    alpha = np.array(alpha_occ)
    beta = np.array(beta_occ)
    tot_orb = np.array(alpha_occ+beta_occ)

    density = orbitals2Density(tot_orb, domain, sigma)
    alpha_density = orbitals2Density(alpha, domain, sigma)
    beta_density = orbitals2Density(beta, domain, sigma)

    return density, alpha_density, beta_density

dom_num = 800
dom = np.linspace(-0.5, 0, dom_num)
sig = 0.02

infile = sys.argv[1]


if os.path.isdir(infile):
    den_list, alpha_list, beta_list = [], [], []
    for filename in os.listdir(infile):
        density, alpha_density, beta_density = genDensity(infile+'/'+filename, dom, sig)
        den_list.append(density)
#        alpha_list.append(alpha_density)
#        beta_list.append(beta_density)
        plt.plot(dom, beta_density, label=filename+'_beta')
        plt.plot(dom, alpha_density, label=filename+'_alpha')
#        plt.plot(dom, density, label=filename)
#        print(filename)

#    plt.plot(dom, den_list[3] - den_list[0], label='nanFe')
#    plt.plot(dom, den_list[1] - den_list[2], label='porFe')
#    den_dif = den_list[3] - den_list[0] - den_list[1] + den_list[2]
#    plt.plot(dom, den_dif, label='porFe')
    

    

    plt.legend()
    plt.show()


elif os.path.isfile(infile):
    density, alpha_density, beta_density = genDensity(infile, dom, sig)
#    plt.plot(dom, density)
    plt.plot(dom, alpha_density, label='alpha')
    plt.plot(dom, beta_density, label='beta')
    plt.legend()
    plt.show()







