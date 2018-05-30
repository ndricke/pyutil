#####################################
#    Lexie McIsaac                  #
#    5/31/16                        #
#    CHEM 268 Final Project         #
#    Hartree-Fock Code              #
#####################################


'''
USAGE:
    
    from the command line, type "python hf.py [molecule] [ngauss] [verbosity] > [output file]"
    if you give no options, it will calculate H2 in STO-3G with verbosity = 0
    you can skip piping it to an output file if you want the output in your terminal 
    
    [molecule] options: 
                        H2    =  the hydrogen molecule
                        HeH   = HeH+ (the helium hydride ion)
                        H3    = the H3+ cation
                        
    [ngauss] options (number of Gaussians in STO basis set):
                        2 = STO-2G
                        3 = STO-3G
                        4 = STO-4G
                        5 = STO-5G
                        6 = STO-6G
                        SO = STO-3G using the helium exponents from Szabo and Ostlund rather than GAMESS
    
    [verbosity] options (optional argument):
                        0    --prints just the electronic energy for each iteration and the final energy (default)
                        1    --also prints the 1-electron integrals and core Hamiltonian
                        2    --also prints the transformation matrix to orthogonalize the basis
                        3    --also prints all the matrices (fock, density, etc) for each SCF iteration
                        4    --also prints all the 2-electron integrals
'''

import numpy as np
import numpy.linalg as npl
import scipy
from scipy import special 
import sys

###############################
#                             # 
#  BASIS SET INFO             #
#                             #
###############################


# first index = gaussian # (not really meaningful, just need to iterate)
# second index: 0 = alpha [ exponent ], 1 = d [ coefficient ]

# from basis set exchange:
H_basis = [[3.42525091,0.15432897],[0.62391373,0.53532814],[0.16885540,0.44463454]]
He_basis = [[6.36242139,0.15432897],[1.15892300,0.53532814],[0.31364979,0.44463454]]

# the sto-3g used by Szabo and Ostlund, different zeta for He, zeta_He = 2.0925:
He_basis_SO = [[9.753934616,0.15432897],[1.776691148,0.53532814],[0.4808442903,0.44463454]]


# put it all together:                                                             
# index = atomic number (just H, He for now)
sto_3g = [H_basis,He_basis]               # standard STO-3G in GAMESS, Gaussian, etc
so_sto_3g = [H_basis, He_basis_SO]        # for comparison with Szabo & Ostlund's sample calculation

# STO-2G
sto2g_h = [[1.3097564,0.430128498301],[0.2331360,0.678913530502]]
sto2g_he = [[2.4328793,0.430128498301],[0.4330513,0.678913530502]]
sto_2g = [sto2g_h, sto2g_he]

# STO-4G
sto4g_h = [[8.0214202,0.056752420805],[1.4678211,0.260141355021],[0.4077768,0.532846114343],
            [0.1353374,0.291625440523]]
sto4g_he = [[14.8998297,0.056752420805],[2.7264853,0.260141355021],[0.7574475,0.532846114343],
            [0.2513900,0.291625440523]]
sto_4g = [sto4g_h, sto4g_he]

# STO-5G
sto5g_h = [[17.3835474,0.022140553120],[3.1854892,0.113541152000],[0.8897299,0.331816148399],
            [0.3037874,0.482570071298],[0.1144785,0.193572196599]]

sto5g_he = [[32.2900297,0.022140553120],[5.9170628,0.113541152000],[1.6526779,0.331816148399],
            [0.5642867,0.482570071298],[0.2126444,0.193572196599]]
sto_5g = [sto5g_h, sto5g_he]

# STO-6G
sto6g_h = [[35.52322122,0.00916359628], [6.513143725,0.04936149294],[1.822142904,0.16853830490],
            [0.625955266,0.37056279970], [0.243076747,0.41649152980], [0.100112428,0.13033408410]]            
sto6g_he = [[65.98456824,0.00916359628],[12.09819836,0.04936149294],[3.384639924,0.16853830490],
            [1.162715163,0.37056279970],[0.451516322,0.41649152980],[0.185959356,0.13033408410]]
sto_6g = [sto6g_h, sto6g_he]

STO = [sto_2g, sto_3g, sto_4g, sto_5g, sto_6g, so_sto_3g] # list of all the basis sets

#################################
#                               #
#      CONSTANT PARAMETERS:     #
#                               #
#################################

BtA = 0.52917721067    # convert from bohr to angstroms
AtB = 1/BtA            # convert A to Bohr

# ATOMIC NUMBERS:
ZHH = [1,1]      # list of atomic numbers for H2 
ZHeH = [1,2]     # list of atomic numbers for HeH+
ZH3 = [1,1,1]    # list of atomic numbers for H3+

# BOND LENGTHS:
R_HH = 1.4       # bond length for H2 (in Bohr), from Szabo & Ostlund
R_HeH = 1.4632   # bond length for HeH+, from Szabo & Ostlund

# CARTESIAN COORDINATES: (in Bohr)
HH_coord = [np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 0.0, R_HH])]
            
HeH_coord = [np.array([0.0, 0.0, 0.0]),
             np.array([0.0, 0.0, R_HeH])]
             
H3_coord = [np.array([0.0, 0.0, 0.0]), 
            np.array([1*AtB, 0.0, 0.0]),
            np.array([.5*AtB,.865*AtB,0.0])] # from DeKock & Gray, Chemical Structure and Bonding, p 272


# DICTIONARIES CONTAINING ALL MOLECULE INFO
H2 = {"name": "H2","atoms": ["H", "H"], "coordinates": HH_coord, "Z": ZHH}
HeH = {"name": "HeH+","atoms": ["H", "He"], "coordinates": HeH_coord, "Z": ZHeH}
H3 = {"name": "H3+","atoms": ["H", "H", "H"], "coordinates": H3_coord, "Z": ZH3}

molecules = {"H2":H2, "HeH":HeH, "H3":H3}



###############################
#                             #
#     PARSE INPUT             #
#                             #
###############################

# choose the molecule:
mol_inp = sys.argv[1]      

if mol_inp in molecules:
    molecule = molecules[mol_inp]
else: # default in case input is bad
    print "ERROR: molecule input specified incorrectly, performing calculation on H2"
    molecule = H2

# choose the basis set
ng = sys.argv[2]
if ng != "SO" and int(ng) < 7 and int(ng) > 1:
    basis = STO[int(sys.argv[2])-2]    
elif ng == "SO":
    basis = STO[-1]
else:  # default in case input is bad
    print "ERROR: basis set specified incorrectly, performing calculation using STO-3G"
    basis = STO[1]


# choose the verbosity
v = sys.argv[3]    
if int(v) < 5 and int(v) > -1:
    verbose = int(v)
elif not v:      # default if not specified
    verbosity = 0
else:            # default if input is bad
    print "ERROR: verbosity specified incorrectly, performing calculation using verbosity = 0"
     
nbasis = len(molecule["atoms"])          # total number of basis functions 
ngauss = len(basis[0])                   # number of gaussians used to make up the CGF basis functions 



######################
#                    #
#  HELPER FUNCTIONS  #
#                    #
######################

def enuc(pos,Z): 
    '''
    Calculate the nuclear repulsion energy
    
    From Szabo & Ostlund p. 150, eq. 3.185
    
    Inputs:
            pos = list of positions of all atoms
            Z = list of atomic numbers for all atoms
    Returns:
            nre = nuclear repulsion energy
    '''
    nre = 0
    natom = len(Z)
    for i in range(0,natom):
        for j in range(i+1,natom):
            R = npl.norm(pos[i]-pos[j]) # bond length
            nre += Z[i]*Z[j]/R
    return nre
    

def integrals(basis_coeff, pos,Zs, nbasis, ngauss): 
    '''
    Calculate the one electron integrals
    
    From Szabo & Ostlund, mostly p. 153-161 and p. 411-415
    
    Inputs:
            basis_coeff = list of coefficients for your basis set, in the form [alpha,d] (for more details on the format, see below in parameters)
            pos = list of positions of the atoms
            Zs = list of atomic numbers of the atoms
            nbasis = total number of basis functions
            ngauss = number of gaussian functions making up each basis function (N in STO-NG)
    
    Returns:
            S = overlap matrix
            T = kinetic energy matrix
            V = list of nucleus-electron potential energy matrices
            (all matrices are numpy arrays)   
    '''
    natom = len(pos)
    
            
    # initialize the matrices
    S = np.zeros((nbasis,nbasis))
    T = np.zeros((nbasis,nbasis))
        
    V = []  # V is potential caused by each atom, so it's a list of initialized matrices
    for i in range(0,natom):
        V.append(np.zeros((nbasis,nbasis)))
    
    
    # mu and nu iterate over the STO basis functions (Contracted Gaussian Fucntions/CGFs)
    for mu in range(0,nbasis):    
        R_mu = pos[mu]      # center of the CGF
        Zmu = Zs[mu]        # atomic number of the atom the CGF describes

        for nu in range(0, nbasis):
            R_nu = pos[nu]  # center of the CGF
            Znu = Zs[nu]    # atomic number of the atom the CGF describes

            Rmn = npl.norm(R_mu-R_nu)     # distance between the centers
            
            # initialize the matrix elements; m = mu, n = nu
            Smn = 0
            Tmn = 0
            
            # p and q iterate over the N Gaussians that make up the STO-NG function
            for p in range(0,ngauss):
                alpha_p = basis_coeff[Zmu-1][p][0]        # exponent 
                norm_p = (2*alpha_p/np.pi)**0.75          # normalization constant
                d_p = basis_coeff[Zmu-1][p][1]*norm_p     # normalized expansion coefficient
                
                for q in range(0,ngauss):
                    alpha_q = basis_coeff[Znu-1][q][0]         # exponent
                    norm_q = (2*alpha_q/np.pi)**0.75           # normalization constant
                    d_q = basis_coeff[Znu-1][q][1]*norm_q      # normalized expansion coefficient

                    # calculate the overlap integral for the p,q pair:
                    Spq = overlap(alpha_p, alpha_q, Rmn)
                    # normalize/multiply by contraction coefficients and add it to the total overlap for the mu,nu pair:
                    Smn += d_q*d_p*Spq

                    # calculate the kinetic energy integral for the p,q pair:
                    Tpq = kinetic(alpha_p, alpha_q, Rmn, Spq)
                    # normalize/multiply by contraction coefficients add it to the total kinetic for the mu,nu pair
                    Tmn += d_q*d_p*Tpq

                    # center of new GF = GF(mu)*GF(nu)    in Szabo & Ostlund this is Rp
                    Rf = (alpha_p*R_mu + alpha_q*R_nu)/(alpha_p+alpha_q)
                    
                    
                    # loop through all the atoms to calculate the potential caused by each one
                    for i in range(0, natom):
                        Rif = npl.norm(Rf - pos[i])       # distance between the new gaussian, F, and the atom
                        Vipq = potential(alpha_p, alpha_q, Rmn, Zs[i],Rif)     # calculate potential energy integral
                        V[i][mu][nu] += d_p*d_q*Vipq      # normalize, etc, and add to matrix
                        
                        

            S[mu][nu] = Smn
            T[mu][nu] = Tmn
             
    return S, T, V



def overlap(alpha_p, alpha_q, R):
    '''
    Calculate the overlap integral between two Gaussians
    
    From Szabo and Ostlund, p 412, eqn A.9
    
    Inputs:
            alpha_p, alpha_q = gaussian exponents for the 2 gaussian functions
            R = distance between the two gaussian centers
    Returns:
            Spq = overlap integral (float)
    '''
    Spq = ((np.pi/(alpha_p + alpha_q))**1.5)*np.exp(-alpha_p*alpha_q*(R**2)/(alpha_p+alpha_q))
    return Spq



def kinetic(alpha_p, alpha_q, R, Spq):
    '''
    Calculate the kinetic energy integral between two Gaussians
    
    Szabo and Ostlund p. 412, eqn A.11
    
    Inputs:
            alpha_p, alpha_q = gaussian exponents for the 2 gaussian functions
            R = distance between the two gaussian centers
            Spq = overlap integral between the two gaussians
    Returns:
            Tpq = kinetic energy integral (float)
    '''        
    Tpq = Spq*(alpha_p*alpha_q/(alpha_p + alpha_q))*(3 - 2*alpha_p*alpha_q*(R**2)/(alpha_p+ alpha_q))
    return Tpq
    
    
    
def potential(alpha_p, alpha_q, Rpq, Zc, Rfc):
    '''
    Calculates the potential energy integrals due to atom C
    
    Szabo and Ostlund p. 415, eqn A.33
    
    Inputs:
            alpha_p, alpha_q = gaussian exponents for the 2 gaussian functions
            Rpq = distance between centers of p and q
            Zc = atomic number of the nucleus causing the potential energy term
            Rfc = distance between the combined gaussian F and atom C
    Returns:
            Vpq = potential energy integral (float)
    '''
    t = (alpha_p+alpha_q)*(Rfc**2) # argument for boys function
    Vpq = (-2*np.pi*Zc/(alpha_p+alpha_q))*np.exp(-alpha_p*alpha_q*(Rpq**2)/(alpha_p+alpha_q))*F0(t)
    return Vpq
      
      
      
def F0(t): 
    '''
    Calculates Boys' function for n = 0
    
    See Szabo and Ostlund p 415 eqn A.32 and Molecular Electronic Structure Theory 
    by Helgaker, Jorgensen, and Olsen eqn 9.8.6
    '''
    if t ==0 : 
        return 1
    else: 
        return 0.5*np.sqrt(np.pi/t)*scipy.special.erf(np.sqrt(t))   
      
           
                     
def get2e(basis_coeff, pos,Zs, ngauss, nbasis):
    '''
    Calculate the 2 electron integrals
    
    See Szabo and Ostlund p 141 and p 416 eqn A.41 for more info
    
    Inputs:
            basis_coeff = list of coefficients for your basis set, in the form [alpha,d] 
                            (for more details on the format, see below in parameters)
            pos = list of position vectors for all atoms
            Zs = list of atomic numbers of all atoms
            ngauss = number of gaussians that make up each basis function
            nbasis = total number of basis functions
    Returns:
            ints = a dictionary. the keys are a tuple of the electron indices (mu,nu,lambda,sigma) 
                    and the entries are floating point numbers of the integral value
    '''

    ints = {}                     # initialize the dictionary
    
    # mu, nu iterate over the STO basis functions
    for mu in range(0,nbasis):    
        Rmu = pos[mu]             # center of the basis function
        Zmu = Zs[mu]              # atomic number of the atom the basis function belongs to
        for nu in range(0,nbasis):
            Rnu = pos[nu]         # center of the basis function
            Znu = Zs[nu]          # atomic number of the atom the basis function belongs to
            
            Rmn = npl.norm(Rmu-Rnu)   # distance between the centers of the two functions
            
            # lambd, sigma iterate over the STO functions in the same way as mu, nu
            for lambd in range(0,nbasis):
                Rl = pos[lambd]
                Zl = Zs[lambd]
                
                for sigma in range(0,nbasis):
                    Rs = pos[sigma]
                    Zsig = Zs[sigma]
                    
                    Rls = npl.norm(Rl - Rs)
                    
                    Vmnls = 0       # initialize the integral value
                    
                    # p, q, r, s iterate over the Gaussians that make up the STO functions
                    for p in range(0,ngauss):
                        alpha_p = basis_coeff[Zmu-1][p][0]        # exponent 
                        norm_p = (2*alpha_p/np.pi)**0.75          # normalization constant
                        d_p = basis_coeff[Zmu-1][p][1]*norm_p     # normalized expansion constant
                        
                        for q in range(0,ngauss):
                            alpha_q = basis_coeff[Znu-1][q][0]        # exponent 
                            norm_q = (2*alpha_q/np.pi)**0.75          # normalization constant
                            d_q = basis_coeff[Znu-1][q][1]*norm_q     # normalized expansion constant                       
                            
                            # center of the new gaussian formed by multiplying p*q (in Szabo and Ostlund this is Rp)
                            Rf = (alpha_p*Rmu + alpha_q*Rnu)/(alpha_p+alpha_q)
                            
                            for r in range(0,ngauss):
                                alpha_r = basis_coeff[Zl-1][r][0]        # exponent 
                                norm_r = (2*alpha_r/np.pi)**0.75         # normalization constant
                                d_r = basis_coeff[Zl-1][r][1]*norm_r     # normalized expansion constant
                                                        
                                for s in range(0,ngauss):
                                    alpha_s = basis_coeff[Zsig-1][s][0]        # exponent 
                                    norm_s = (2*alpha_s/np.pi)**0.75           # normalization constant
                                    d_s = basis_coeff[Zsig-1][s][1]*norm_s     # normalized expansion constant                                   
                                    
                                    # center of the new gaussian formed by multiplying r*s (in S&O this is Rq)
                                    Rg = (alpha_r*Rl + alpha_s*Rs)/(alpha_r+alpha_s)
                                    Rfg = npl.norm(Rf - Rg)  # distance between F and G
                                    
                                    # find the potential energy from the pqrs term
                                    Vpqrs = calc2e(alpha_p, alpha_q, alpha_r, alpha_s, Rmn, Rls, Rfg)
                                    Vmnls += d_p*d_q*d_r*d_s*Vpqrs   # normalize and add to the mnls integral value
                    
                    mnls = (mu+1,nu+1,lambd+1,sigma+1)   # create the dictionary key
                    ints[mnls]= Vmnls                    # add the integral to the dictionary
    return ints
                                    
                                    

def calc2e(alpha_p, alpha_q, alpha_r, alpha_s, Rpq, Rrs, Rfg):
    '''
    Calculates the 2-electron integral between four gaussians
    
    See Szabo and Ostlund p 416 eqn A.41
    
    Inputs:
            alpha_p, alpha_q, alpha_r, alpha_s = exponents for the gaussians
            Rpq = distance between the centers of p and q
            Rrs = distance between the centers of r and s
            Rfg = distance between the new gaussians formed when multiplying p*q and r*s
    Returns:
            teint = the 2 electron integral for pqrs, a floating point number
    '''
    pref = 2*(np.pi**2.5)/((alpha_p + alpha_q)*(alpha_r + alpha_s)*np.sqrt(alpha_p + alpha_q + alpha_r + alpha_s))
    exp = np.exp((-alpha_p*alpha_q*(Rpq**2)/(alpha_p + alpha_q)) - (alpha_r*alpha_s*(Rrs**2)/(alpha_r+alpha_s)))
    arg = (alpha_p + alpha_q)*(alpha_r + alpha_s)*(Rfg**2)/(alpha_p + alpha_q + alpha_r + alpha_s)
    teint = pref*exp*F0(arg)
    return teint
  
  
  
def transformation(S):
    '''
    Find the transformation vector, X, that orthogonalizes the basis set
    Uses symmetric orthogonalization, as described in Szabo and Ostlund p 143
    
    Inputs:
        S = overlap matrix
    Returns:
        X = transformation matrix (S^-.5)
        X_dagger = conjugate transpose of X
    '''
    S_eval, S_evec = npl.eigh(S)                     # extract eigenvalues and eigenvectors from overlap matrix
                                                    
    s_sqrt = np.zeros((len(S_eval), len(S_eval)))    # initialize s^-0.5, the diagonalized overlap matrix to the -1/2 power
    for i in range(0,len(S_eval)):                                 # form s^-0.5
        s_sqrt[i][i] = S_eval[i]**(-.5)   
    
    U = S_evec                                   # find the unitary transform matrix
    U_dagger = np.transpose(U)                        # find the conjugate transpose of the unitary matrix
    X = np.dot(U, np.dot(s_sqrt, U_dagger))      # form X = S^-0.5, the transform matrix to orthogonalize the basis set
    X_dagger = np.transpose(X)                   # conjugate transpose is just transpose since all values are real
    return X, X_dagger



def makeP(C, numatoms, numbasis):
    '''
    Calculates the density matrix from the orbital expansion coefficients
    
    See Szabo and Ostlund p 139 eqn 3.145
    
    Inputs:
            C = matrix with orbital expansion coefficients
            numatoms = number of atoms in the system
            numbasis = total number of basis functions
    Returns:
            P = numpy array with the density matrix
    '''
    P = np.zeros((numbasis, numbasis))     # initialize
    
    # mu, nu iterate over STO basis functions
    for mu in range(0, numbasis):
        for nu in range(0, numbasis):
            Pmn = 0
            
            # a iterates over the atoms
            for a in range(0, numatoms/2):
                Pa = C[mu][a]*C[nu][a]
                Pmn += 2*Pa
            P[mu][nu] = Pmn
    return P



def makeG(P, numbasis, twoes):
    '''
    Calculates the G matrix, a helper matrix for constructing the Fock matrix
    
    See Szabo and Ostlund p 141, eqn 3.154
    
    Inputs: 
            P = density matrix (a numpy array)
            numbasis = total number of basis functions
            twoes = dictionary with 2-electron integrals
    Returns:
            G = numpy array with the G matrix
    '''
    G = np.zeros((numbasis,numbasis))            # initialize the matrix
    
    if P.all() == G.all():                       # if P is all 0's, just return the initial all zero G matrix
        return G
    
    for mu in range(0, numbasis):                # mu, nu iterates through the basis functions
        for nu in range(0, numbasis):
            Gmn = 0                              # initialize G[mu][nu]
            for sigma in range(0, numbasis):     # sigma, lambd iterates through the basis functions again (to get 2 e- integrals)
                for lambd in range(0, numbasis):
                    # calculate the contribution to G for the sigma, lambda pair
                    Gls = P[sigma][lambd] * (twoes[(mu+1, nu+1, sigma+1, lambd+1)] - 0.5*twoes[(mu+1, lambd+1, sigma+1, nu+1)])
                    Gmn += Gls   # add it to G[mu][nu]

            G[mu][nu] = Gmn    
    return G



def getEel(P, Hcore, F, nbasis):
    '''
    Calculates the electronic energy
    
    See Szabo and Ostlund p 150 eqn 3.184
    
    Inputs:
            P = density matrix (numpy array)
            Hcore = core hamiltonian matrix (numpy array)
            F = fock matrix (numpy array)
            nbasis = total number of basis functions
    Returns:
            E0 = electronic energy (float)
    '''
    E0 = 0
    for mu in range(0,nbasis):
        for nu in range(0, nbasis):
            E0 += 0.5*P[mu][nu]*(Hcore[mu][nu]+F[mu][nu]) 
    return E0     




################################
#                              #
#    HARTREE-FOCK PROGRAM:     #
#                              #
################################

def main(molecule, nbasis, ngauss):
    '''
    Main function for the Hartree-Fock SCF method
    
    Inputs:
        molecule = dictionary with all molecule info
        nbasis = total number of basis functions
        ngauss = number of gaussians per basis function
    '''
    
    name = molecule["name"]           # name of the molecule
    atoms = molecule["atoms"]         # names of the atoms 
    pos = molecule["coordinates"]     # position vector to use in calculation 
    Z = molecule["Z"]                 # atomic numbers to use in calculation
    natoms = len(atoms)               # number of atoms in the molecule
    
    print "Hartree-Fock calculation for", name, "in STO-{}G".format(ngauss)
    print ""
    print "Molecular Geometry (Bohr):"
    print "{:13} {:7} {:7} {:7}".format("Atom", "X", "Y", "Z")
    for i in range(0,natoms):
        print "{:8} {:7.5} {:7.5} {:7.5}".format(atoms[i], pos[i][0], pos[i][1], pos[i][2])
    print ""
    print ""
        
       
         
    '''
    GET INTEGRALS
    '''       
    # get one electron integrals                                                                             
    S,T,V  = integrals(basis, pos, Z, nbasis, ngauss)    
    
    # sum the potential energies due to each atom into one V term
    Vtot = np.zeros((nbasis,nbasis))
    for Vc in V:
        Vtot += Vc
    
    Hcore = T + Vtot     # form core hamiltonian
    
    # get 2 e- integrals
    twoes = get2e(basis, pos, Z,ngauss,nbasis)       
    
    
    if verbose > 0:
        print "1-ELECTRON INTEGRALS:"
        print "S (Overlap Matrix):"
        print  S
        print ""
        print "T (Kinetic Energy Matrix):"
        print T
        print ""
        print "V (Electron-Nucleus Potential Matrices):"
        for i,Vi in enumerate(V):
            print "Atom #", i+1, ":"
            print Vi
        print ""
        print "Core Hamiltonian:"
        print Hcore
        print ""
    if verbose > 3:
        print "2-ELECTRON INTEGRALS"
        for i in twoes: print i, twoes[i]
        print ""
    
    
    '''
    ORTHOGONALIZE THE BASIS
    '''
    # Get the transformation matrix to orthogonalize the basis set (X = S^-.5)
    X, X_dagger = transformation(S)
    
    if verbose > 1:
        print "X (Transformation Matrix):"
        print X
        print ""
    
    
    '''
    PERFORM SCF ITERATIONS
    
    See Szabo and Ostlund p 145-150
    '''
    print "Beginning SCF Procedure"
    print ""
    
    P = np.zeros((nbasis,nbasis))         # initialize P, the density matrix
    
    maxit = 50                            # maximum number of SCF iterations
    conv = 10**-6                         # convergence criterion for the energy 
    delta = 10                            # set delta (= change in energy) to an arbitrary number larger than conv
    Eel = 0                              # set electronic energy to arbitrary large number
    n = 0                                 # n = # of iterations

    while delta > conv and n <= maxit:
        n += 1                            # keep track of iterations for bookkeeping purposes 
        print "Iteration number ", n
        G = makeG(P, nbasis, twoes)       # calculate the G matrix
        
        F = Hcore + G                     # form the fock matrix, F
        
        Fprime = np.dot(X_dagger,np.dot(F,X))     # F' = X_dagger*F*X is the Fock matrix in the orthogonalized basis
    
        E, Cprime = npl.eigh(Fprime)              # diagonalize F' to get C', the orbital expansion coefficients in the orthogonalized basis
                                                # E is the orbital energies in the orthogonal basis, but isn't really used
        
        C = np.dot(X, Cprime)                     # C = X*C' yields the orbital expansion coefficients
        
        P = makeP(C, natoms, nbasis)              # calculate the new density matrix from the new values of C
        
        Eold = Eel                          # save the electronic energy from the last iteration
        Eel = getEel(P, Hcore, F, nbasis)   # calculate the electronic energy for this iteration
    
        delta = np.abs(Eel - Eold)          # check for convergence
        
        if verbose > 2:
            print "G matrix:"
            print G
            print "F matrix:"
            print F
            print "F' matrix:"
            print Fprime
            print "E matrix:"
            print E
            print "C' matrix:"
            print Cprime
            print "C matrix:"
            print C
            print "Density Matrix:"
            print P
            
        print "Electronic energy = ", Eel, "Hartrees"
        print "Change in energy = ", delta, "Hartrees"
        print ""
    
    print ""
    print "ENERGY CONVERGED AFTER", n, "ITERATIONS"
    print ""
    
    
    Enuc = enuc(pos, Z)
    Etot = Eel + Enuc
    
    print "Final electronic energy = ", Eel, "Hartrees"
    print "Nuclear reupulsion energy =", Enuc, "Hartrees"
    print "Final total energy = ", Etot, "Hartrees"
    print ""
    print ""
    print ""
    print "MOLECULAR ORBITAL INFO:"
    print "Orbital Energies:", E
    print "Expansion Coefficients:"
    print C

        
        
        
if __name__=='__main__':
	main(molecule, nbasis, ngauss)

    









    
