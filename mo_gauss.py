import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

import Qdata

font = {'size':14}
mpl.rc('font',**font)



def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def orbitals2Density(mo_energies, dom, sig):
    density = np.zeros(dom_num)
    for mu in mo_energies:
        #print(mu)
        density += gaussian(dom, mu, sig)
    return density

def genDensity(infile, domain, sigma=0.005):
    qdata = Qdata.Qdata()
    qdata.readFile(infile)

    alpha_occ = []
    for item in qdata.alpha_occ:
        try: alpha_occ.append(float(item))
        except: pass
    beta_occ = []
    for item in qdata.beta_occ:
        try: beta_occ.append(float(item))
        except: pass
    #alpha_occ = [float(item) for item in qdata.alpha_occ]
    #beta_occ = [float(item) for item in qdata.beta_occ]
    alpha = np.array(alpha_occ)
    beta = np.array(beta_occ)
    tot_orb = np.array(alpha_occ+beta_occ)

    #density = orbitals2Density(tot_orb, domain, sigma)
    alpha_density = orbitals2Density(alpha, domain, sigma)
    beta_density = orbitals2Density(beta, domain, sigma)
    density = alpha_density + beta_density

    return density, alpha_density, beta_density

dom_num = 30000
#x_range = [-15.,0.25]
#x_range = [-10.22,-10.1]
#x_range = [-20.,0.]
x_range = [-1.03, -0.12]
#x_range = [-19.25,-19.12]
dom = np.linspace(x_range[0], x_range[1], dom_num)
#sig = 0.05 # 0.2 is a decent depiction of overall band shift
sig = 0.01
#sig = 0.0001
#sig = 0.001

infile = sys.argv[1]


if os.path.isdir(infile):
    den_list, alpha_list, beta_list = [], [], []
    for filename in os.listdir(infile):
        print(filename)
        density, alpha_density, beta_density = genDensity(infile+'/'+filename, dom, sig)
        den_list.append(density)
#        alpha_list.append(alpha_density)
#        beta_list.append(beta_density)
#        plt.plot(dom, beta_density, label=filename+'_beta')
#        plt.plot(dom, alpha_density, label=filename+'_alpha')
#        plt.plot(dom, density, label=filename)
#        print(filename)

    #print(den_list[1])

    #plt.plot(dom, den_list[4] - den_list[0], label='nanFe')
    #plt.plot(dom, den_list[2] - den_list[3], label='porFe')
    #plt.plot(dom, den_list[1], label='O2')
    plt.plot(dom, den_list[4] - den_list[0] - den_list[1], label='nanFe')
    plt.plot(dom, den_list[2] - den_list[3] - den_list[1], label='porFe')

    #den_dif = den_list[4] - den_list[0] - den_list[2] + den_list[3]
    #plt.plot(dom, den_dif, label='nanFeDiff - porFeDiff')




    plt.plot(x_range,[0.,0.], color='k')
    plt.xlim(x_range)
    plt.xlabel("Band Energy (Ht)")
    plt.ylabel("Density Difference")
    plt.legend()
    plt.show()
    #plt.savefig("nanporFe_totdiff_s0p01.png", transparent=True, bbox_inches='tight', pad_inches=0.05)


elif os.path.isfile(infile):
    density, alpha_density, beta_density = genDensity(infile, dom, sig)
    plt.plot(dom, density, label = 'Total')
    #plt.plot(dom, alpha_density, label='alpha')
    #plt.plot(dom, beta_density, label='beta')
    
    plt.xlabel("Band Energy (Ht)")
    plt.ylabel("Density Difference")
    plt.legend()
    plt.show()




