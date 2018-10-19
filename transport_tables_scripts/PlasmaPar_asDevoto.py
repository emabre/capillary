import numpy as np
import constantsGAU_ema as gau
import ionization as iz
import itertools as it

# <codecell>

# alpha, beta, A, B : [l,s,m,p]
A = np.zeros((3,5,3,3))
B = np.zeros((1,5,3,3))

# For 6th approx
#A = np.zeros((4,11,6,6))
#B = np.zeros((1,11,6,6))

# q00
B[0,0,0,0] = 8*1
# q01
B[0,0,0,1] = 8*5/2; B[0,1,0,1] = 8*(-3)
# q11
A[1,1,1,1] = 8*1
B[0,0,1,1] = 8*25/4; B[0,1,1,1] = 8*(-15); B[0,2,1,1] = 8*12
# q02
B[0,0,0,2] = 8*35/8; B[0,1,0,2] = 8*(-21/2); B[0,2,0,2] = 8*6
# q12
A[1,1,1,2] = 8*7/4; A[1,2,1,2] = 8*(-2)
B[0,0,1,2] = 8*175/16; B[0,1,1,2] = 8*-315/8; B[0,2,1,2] = 8*57;
B[0,3,1,2] = 8*-30
# q22
A[1,1,2,2] = 8*77/16; A[1,2,2,2] = 8*-7; A[1,3,2,2] = 8*5
B[0,0,2,2] = 8*1225/64; B[0,1,2,2] = 8*-735/8; B[0,2,2,2] = 8*399/2;
B[0,3,2,2] = 8*-210; B[0,4,2,2] = 8*90
## q03
#B[0,0,0,3] = 8*105/16; B[0,1,0,3] = 8*-189/8; B[0,2,0,3] = 8*27;
#B[0,3,0,3] = 8*-10
## q13
#A[1,1,1,3] = 8*63/32; A[1,2,1,3] = 8*-9/2; A[1,3,1,3] = 8*5/2
#B[0,0,1,3] = -8*525/32; B[0,1,1,3] = -8*-315/4; B[0,2,1,3] = -8*162;
#B[0,3,1,3] = -8*-160; B[0,4,1,3] =-8*60 
## q23
#A[1,1,2,3] = 8*945/128; A[1,2,2,3] = 8*-261/16; A[1,3,2,3] = 8*125/8;
#A[1,4,2,3] = 8*-15/2
#B[0,0,2,3] = 8*3675/128; B[0,1,2,3] = 8*-11025/64; B[0,2,2,3] = 8*1953/4
#B[0,3,2,3] = 8*-1505/2; B[0,4,2,3] = 8*615; B[0,5,2,3] = 8*-210
## q33
#A[1,1,3,3] = 8*14553/1024; A[1,2,3,3] = 8*-1215/32; A[1,3,3,3] = 8*1565/32;
#A[1,4,3,3] = 8*-135/4; A[1,5,3,3] = 8*105/8; A[3,3,3,3] = 8*1;
#B[0,0,3,3] = 8*11025/256; B[0,1,2,3] = 8*-19845/64; B[0,2,2,3] = 8*17577/16;
#B[0,3,2,3] = 8*-4515/2; B[0,4,2,3] = 8*5535/2; B[0,5,2,3] = 8*-1890
#B[0,6,2,3] = 8*560

# <codecell>
# g_i parameters, from table VI of Bruno, Phys Plasmas 17, 112315 (2010)
g_eH = {}
g_eH[(1,1)] = [10.35291134,   -1.58301162, 12.45844036, -0.23285190,    5.36628573e-2, -5.34372929, 9.35561752,  -2.15463427]
g_eH[(1,2)] = [10.09959930,   -1.50068352, 12.54524872, -7.29529868e-2, 4.37301268e-2, -5.73166847, 9.09798179,  -2.13265127]
g_eH[(1,3)] = [9.84443341,    -1.42568830, 12.97194554, -9.24067489e-2, 2.32754987e-2, -5.71948057, 8.8325970,   -2.05797013]
g_eH[(1,4)] = [9.65959253,    -1.38293428, 13.18865311, -6.34310879e-2, 1.11968653e-2, -5.77977010, 8.64525855,  -2.02715634]
g_eH[(1,5)] = [9.50113220,    -1.36043711, 13.23885240, -4.44172748e-2, 7.88864536e-3, -5.83593794, 8.49000000,  -2.01418763]
g_eH[(2,2)] = [10.33445440,   -1.44880911, 12.08534341, -1.86163984e-2, 3.30723152e-2, -6.45649277, 9.15932646,  -2.13494419]
g_eH[(2,3)] = [10.08612484,   -1.39070408, 12.39823810, -3.26532858e-2, 1.75287406e-2, -6.48186913, 8.90783546,  -2.08119318]
g_eH[(2,4)] = [9.89312188,    -1.34820033, 12.63138071, -1.96030794e-2, 4.55766915e-3, -6.50687636, 8.71405704,  -2.04690115]
g_eH[(3,3)] = [9.99023294,    -1.41457896, 12.53875975, -2.82169003e-2, 2.70547916e-2, -6.11376507, 8.89657156,  -2.09400530]

# a_i parameters, from table I of Bruno, Phys Plasmas 17, 112315 (2010)
a_HH = {}
a_HH[(1,1)] = [15.09506044, -1.25710008, 9.57839369, -3.80371463, 0.98646613, 9.25705877, -0.93611707]
a_HH[(1,2)] = [14.14566908, -1.17057105, 9.02830724, -3.00779776, 0.74653903, 9.10299040, -0.68184353]
a_HH[(1,3)] = [13.39722075, -1.09886403, 8.50097335, -2.86025395, 0.85345727, 8.90666490, -0.67571329]
a_HH[(1,4)] = [12.97073246, -1.06479185, 8.18885522, -2.78105132, 0.89401865, 8.73403138, -0.65658782]
a_HH[(1,5)] = [12.69248000, -1.04857945, 7.97861283, -2.73621289, 0.90816787, 8.57840253, -0.63732002]
a_HH[(2,2)] = [22.08948804, -1.85066626, 8.50932055, -7.66943974, 0.77454531, 9.69545318, -0.62104466]
a_HH[(2,3)] = [17.94703897, -1.42488999, 7.66669340, -4.76239721, 1.26783524, 9.53716768, -0.73914215]
a_HH[(2,4)] = [18.78590499, -1.59291967, 7.97734302, -5.66814860, 1.01816360, 9.32328437, -0.60882006]
a_HH[(3,3)] = [13.82986524, -1.01454290, 7.48970759, -3.27628187, 2.08225623, 9.21388055, -1.32086596]

# a_i parameters, from table III of Bruno, Phys Plasmas 17, 112315 (2010)
a_HpH = {}
a_HpH[(1,1)] = [46.68783791, -0.33303803, 4.25686770, -2.03851201, 14.98170958, 8.59618369, -1.65616736]
a_HpH[(1,2)] = [46.68783791, -0.33303803, 3.92217635, -2.00886829, 14.98170958, 8.24501842, -1.65616736]
a_HpH[(1,3)] = [46.68783791, -0.33303803, 3.65740159, -2.01434735, 14.98170958, 7.97534885, -1.65616736]
a_HpH[(1,4)] = [46.68783791, -0.33303803, 3.43102576, -2.04002032, 14.98170958, 7.76086951, -1.65616736]
a_HpH[(1,5)] = [46.68783791, -0.33303803, 3.23079831, -2.07543755, 14.98170958, 7.58613195, -1.65616736]
a_HpH[(2,2)] = [46.68783791, -0.33303803, 4.10212490, -1.85454858, 14.98170958, 8.86285119, -1.65616736]
a_HpH[(2,3)] = [46.68783791, -0.33303803, 3.89701552, -1.76267951, 14.98170958, 8.61913831, -1.65616736]
a_HpH[(2,4)] = [46.68783791, -0.33303803, 3.73496748, -1.69596577, 14.98170958, 8.41234103, -1.65616736]
a_HpH[(3,3)] = [46.68783791, -0.33303803, 3.95678840, -2.00381603, 15.72150840, 8.30656354, -1.79347178]

# d_i parameters, from table V of Bruno, Phys Plasmas 17, 112315 (2010)
d_HHp = {}
d_HHp[(1,1)] = [ 63.5437, -5.0093, 9.8797e-2]
d_HHp[(1,2)] = [ 61.8730, -4.9431, 9.8766e-2]
d_HHp[(1,3)] = [ 60.6364, -4.8936, 9.8767e-2]
d_HHp[(1,4)] = [ 59.6591, -4.8544, 9.8785e-2]
d_HHp[(1,5)] = [ 58.8493, -4.8213, 9.8776e-2]

# <codecell>

def elRes_norm(z,rho,kT, ne=0):
    ''' Electrical resistivity of hydrogen as computed by Devoto in "Simplified expressions
        for the transport properties of ionized monatomic gases"(1967).
        See eq. 16, and take the 3rd approximation for D11.
    '''
    # This should be fixed later!
    if ne==0:
        ne = iz.elec_dens(rho,z)
    else:  # Just a reminder
        print("Fix this crap!")
    
    n_p = rho/(gau.me+gau.mp)
    n = np.array([ne, ne, n_p-ne])
    Q = makeQ(kT, ne)
    # Total number density of free particles in gas (ne + nH+ + nH), assuming 
    # complete dissociation:
    # Ordinary diff coefficient, 3rd approx.
    D = diffu_ee(kT, rho, n, Q)
    eta = rho * kT / ( gau.qe**2 * ne * n.sum() * gau.me * D)
    return eta

def elRes_norm_2(z,rho,kT):
    ''' El. resistivity of hydrogen as in paper "Transport coeff...",Devoto(1966),
        excluding ion current contribution.
        See eq. (29)
    '''
    ne = iz.elec_dens(rho,z)
    Q = makeQ(kT, ne)
    # Total number density of free particles in gas (ne + nH+ + nH), assuming 
    # complete dissociation:
    n_p = rho/(gau.me+gau.mp)
    n = np.array([ne, ne, n_p-ne])
    m = np.array([gau.me, gau.mp, gau.mp+gau.me])
    
    D = diffu(kT, rho, n, m, Q)
    Zi = np.array([-1.0,1.0,0.0])
    sigma = gau.qe**2*n.sum()/(rho*kT) * np.sum(n[1:]*m[1:]*Zi[1:]*D[0,1:])
    return 1/sigma

def elRes_norm_3(z,rho,kT):
    ''' Corrected version of el. resistivity of hydrogen as in paper "Transport coeff...",Devoto(1966),
        including ion current contribution.
        See eq. (28)
    '''
    ne = iz.elec_dens(rho,z)
    Q = makeQ(kT, ne)
    # Total number density of free particles in gas (ne + nH+ + nH), assuming 
    # complete dissociation:
    n_p = rho/(gau.me+gau.mp)
    n = np.array([ne, ne, n_p-ne])
    m = np.array([gau.me, gau.mp, gau.mp+gau.me])
    
    D = diffu(kT, rho, n, m, Q)
    Zi = np.array([-1.0, 1.0, 0.0])
    
    sigma = n[1]*m[1]*Zi[1]*D[0,1] - Zi[1]*n[0]*m[0]*Zi[0]*D[1,0]
    sigma *= gau.qe**2*n.sum()/(rho*kT)
    
    return 1/sigma

def thermCond_e_norm(z, rho, kT, corr=True, B=0.0):
    ''' Translational thermal conductivity of electrons, as computed by Devoto in "Simplified expressions
        for the transport properties of ionized monatomic gases"(1967).
        See eq. 21.
    '''
    ne = iz.elec_dens(rho,z)
    Q = makeQ(kT, ne)
    # Total number density of free particles in gas (ne + nH+ + nH), assuming 
    # complete dissociation:
    n_p = rho/(gau.me+gau.mp)
    n = np.array([ne, ne, n_p-ne])
    
    q11 = q_mp_simple(1, 1, n, Q)
    q12 = q_mp_simple(1, 2, n, Q)
    q22 = q_mp_simple(2, 2, n, Q)
    
    k = 75*ne**2*gau.kB/8 * np.sqrt(2*np.pi*kT/gau.me) * (q11-(q12**2)/q22)**-1
    
    return k

def thermCond_h_norm(z, rho, kT, corr=True, B=0.0):
    ''' Translational thermal conductivity (for H,H+ mixture gas (equilibrium hydrogen,
        but forgetting the free electrons)) as in Devoto,
        "Transport properties of ionized monatomic gases"(1966) (see eq. (17))
    '''
    ne = iz.elec_dens(rho,z)
    # I cannot keep all the collision integrals, as I have no e-
    Q = makeQ(kT, ne)[:,:,1:,1:]
    # Total number density of free particles in gas (ne + nH+ + nH), assuming 
    # complete dissociation:
    n_p = rho/(gau.me+gau.mp)
    # number densities of gas components (H, H+)
    n = np.array([ne, n_p-ne])
    # masses of gas components (H, H+)
    m = np.array([gau.mp, gau.me+gau.mp])
    
    lprime = lambdaPrime(kT, n, m, Q)
    
    D = diffu(kT, rho, n, m, Q)
    E = Eij(D, m)
    DT = thermDiffu(kT, n, m, Q)
    
    addend = 0
    idx = np.arange(len(n))
    for ii,jj in it.product(idx,idx):
        addend += E[ii,jj]*DT[ii]*DT[jj]/(n[ii]*m[ii]*m[jj])
    addend *= rho*gau.kB/n.sum()
    
    return addend+lprime

def thermCond_norm(z, rho, kT):
    ''' Traslational thermal conductivity (for hydrogen gas) as in Devoto, "Transport properties of ionized
        monatomic gases"(1966) (see eq. (17))
    '''
    ne = iz.elec_dens(rho,z)
    Q = makeQ(kT, ne)
    # Total number density of free particles in gas (ne + nH+ + nH), assuming 
    # complete dissociation:
    n_p = rho/(gau.me+gau.mp)
    n = np.array([ne, ne, n_p-ne])
    m = np.array([gau.me, gau.mp, gau.me+gau.mp])
    
    lprime = lambdaPrime(kT, n, m, Q)
    
    D = diffu(kT, rho, n, m, Q)
    E = Eij(D, m)
    DT = thermDiffu(kT, n, m, Q)
    
    addend = 0
    idx = np.arange(len(n))
    for ii,jj in it.product(idx,idx):
        addend += E[ii,jj]*DT[ii]*DT[jj]/(n[ii]*m[ii]*m[jj])
    addend *= rho*gau.kB/n.sum()
    
    return addend+lprime

def thermCond_r_norm(z, rho, kT):
    ''' Reactive thermal conductivity (for hydrogen gas) as in Devoto, "Transport properties of ionized
        monatomic gases"(1966) (see eq. (17))
    '''
    T = kT/gau.kB
    ne = iz.elec_dens(rho,z)
    Q = makeQ(kT, ne)
    # Total number density of free particles in gas (ne + nH+ + nH), assuming 
    # complete dissociation:
    n_p = rho/(gau.me+gau.mp)
    n = np.array([ne, ne, n_p-ne])
    m = np.array([gau.me, gau.mp, gau.me+gau.mp])
    
    # THIS IS PROBABLY WRONG, MAYBE IT SHOULD BE A BIT HIGHER!
    dh = 13.6*1.6e-12  # Entalpy difference in erg
    
    DA, DAT = diffuAmbANDthermDiffuAmb(kT, rho, n, m, Q)
    
    lambda_r = dh * (n.sum()*gau.mp/(2*rho*kT*T)*(n[0]*n[2])/(n[0]+n[2])*dh*DA \
                     + DAT/(gau.mp+T))
    
    return lambda_r

def diffu_ee(kT, rho, n, Q):
    ''' Ordinary diffusion coefficient for electrons,
        3rd approx as in Devoto, "Simplified expressions
        for the transport properties of ionized monatomic gases"(1967).
        n: must contain species number density, with:
            n[0]:electron density
            n[1]:H+ density
            n[2]:H density
    '''
    q00 = q_mp_simple(0, 0, n, Q)
    q01 = q_mp_simple(0, 1, n, Q)
    q11 = q_mp_simple(1, 1, n, Q)
    q12 = q_mp_simple(1, 2, n, Q)
    q22 = q_mp_simple(2, 2, n, Q)
    q02 = q_mp_simple(0, 2, n, Q)
    
#    det_num = q11*q22 - q12**2
#    det_denom = q00 * (q11*q22 - q12**2) + \
#                q01 * (q12*q02 - q01*q22) + \
#                q02 * (q01*q12 - q11*q02)

    num = np.array([[q11, q12],
                    [q12, q22]])
    denom = np.array([[q00, q01, q02],
                      [q01, q11, q12],
                      [q02, q12, q22]])
    
#    D11 = 3*n[0]*rho/(2*n.sum()*gau.me) * np.sqrt(2*np.pi*kT/gau.me)*det_num/det_denom
    D11 = 3*n[0]*rho/(2*n.sum()*gau.me) * np.sqrt(2*np.pi*kT/gau.me)*np.linalg.det(num)/np.linalg.det(denom)
    return D11    

def diffu(kT, rho, n, m, Q):
    ''' Ordinary diffusion coefficients matrix (general, for all the species),
        3rd approx as in Devoto, "Transport properties of Ionized
        monatomic gases"(1966).
        n: array with species number density
        m: array with scecies atomic masses (in same order as n)
    '''
    Nspecs = len(n)
    
    denom = denomD_3(n, m, Q)
    
    num = np.zeros((denom.shape[0]+1, denom.shape[1]+1))
    num[0:denom.shape[0],
        0:denom.shape[1]] = denom.copy()
        
    idx = np.arange(Nspecs)
    D = np.zeros((Nspecs,Nspecs))
    for ii,jj in it.product(idx,idx):
        num[-1,0:Nspecs] = (idx==ii).astype(np.int64)
        num[0:Nspecs,-1] = (idx==jj).astype(np.int64) - (idx==ii).astype(np.int64)
#        print("num:")
#        print(num)
        D[ii,jj] = 3*rho*n[ii]/(2*n.sum()*m[jj])*np.sqrt(2*np.pi*kT/m[ii]) \
            * np.linalg.det(num) / np.linalg.det(denom)
    return D

def thermDiffu(kT, n, m, Q):
    ''' Thermal diffusion coefficients matrix (general, for all the species),
        3rd approx as in Devoto, "Transport properties of Ionized
        n: array with species number density
        m: array with scecies atomic masses (in same order as n)
    '''
    Nspecs = len(n)    
    denom = denomD_3(n, m, Q)
    
    num = np.zeros((denom.shape[0]+1, denom.shape[1]+1))
    num[0:denom.shape[0],
        0:denom.shape[1]] = denom.copy()
        
    idx = np.arange(Nspecs)
    DT = np.zeros(Nspecs)
    for ii in idx:
        num[-1,0:Nspecs] = (idx==ii).astype(np.int64)
        num[Nspecs:2*Nspecs,-1] = n
        DT[ii] = 15*n[ii]*np.sqrt(2*np.pi*m[ii]*kT)/4 \
                 * np.linalg.det(num) / np.linalg.det(denom)
    return DT

def lambdaPrime(kT, n, m, Q):
    ''' Coefficient λ' (see eq 14, "Transport properties...",Devoto(1966))
        computed to the third approx.
    '''
    Nspecs = len(n)
    denom = denomD_3(n, m, Q)
    
    num = np.zeros((denom.shape[0]+1, denom.shape[1]+1))
    num[-1,Nspecs:2*Nspecs] = n/np.sqrt(m)
    num[Nspecs:2*Nspecs,-1] = n
    num[0:denom.shape[0],
        0:denom.shape[1]] = denom.copy()
    
    lprime = -75*gau.kB/8 * np.sqrt(2*np.pi*kT) \
             * np.linalg.det(num) / np.linalg.det(denom)
    return lprime

def diffuAmbANDthermDiffuAmb(kT, rho, n, m, Q):
    ''' Ordinary and thermal ambipolar diffusion coefficients for a partially ionized
        monatomic gas (1 ion species, 1 neutral species, and electrons. (e.g. hydrogen));
        as in Devoto, Phys Fluids 9, 1230 (1966) (see eqs.43,44)
        n: species number densities, with:
            n[0] : electrons
            n[1] : ions
            n[2] : neutral atoms
        analogous is for m and for the shape of Q
        It is important that in place 0 there are electrons, as the simplification
        of terms of the order of m[0]/m[1] is carried out! (see text after eq. (44))
    '''
    if not(3==len(n)==len(m)==Q.shape[2]==Q.shape[3]):
        raise ValueError("Some dimensions mismatch! You should have electrons, \
                         1 ion species, 1 neutral species (in this order!) (e.g. hydrogen)")
    
    D = diffu(kT, rho, n, m, Q)
    DT = thermDiffu(kT, n, m, Q)
    
    DA = 2*D[1,2]*(1 - (m[0]*D[1,0])/(m[1]*D[1,2])*(1-D[0,2]/D[0,1]) )
    DAT = DT[1]*(1 + D[1,0]/D[0,1]*(DT[0]/DT[1] - m[0]/m[1]))
    
    return DA, DAT

def Eij(D, m):
    ''' Matrix E_ij, whose elements appear in equation (17) of papar by Devoto:
        "Transport coefficients..."(1966)
    '''
    E = D.copy()
    for jj in range(D.shape[1]):
        E[:,jj] = D[:,jj]*m[jj]
    E = np.linalg.inv(E)
    return E

def denomD_3(n, m, Q):
    ''' Build matrix q used Devoto (3rd approx), "Transport properties of Ionized
        monatomic gases"(1966) eqs: (8),(9),(14)
        n: array with species number density
        m: array with scecies atomic masses (in same order as n)
    '''
    Nspecs = len(n)
    denom = np.zeros((Nspecs*3,Nspecs*3))
    for mm,pp in it.product([0,1,2],[0,1,2]):
        denom[mm*Nspecs:(mm+1)*Nspecs,
              pp*Nspecs:(pp+1)*Nspecs] = q_mp_complete(mm, pp, n, m, Q)
    return denom

def D_Devoto2Bruno(n,m,D):
    ''' Function to convert diffusion parameters from Devoto's convention
        to Bruno's (Phys.Plasmas14,022303(2007)) convention.
        Note: Bruno's D^k_i is D[i,k]
    '''
    rho = np.sum(n*m)
    n_tot = n.sum()
    D_Bruno = np.empty_like(D)
    for ii in range(D.shape[0]):
        D_Bruno[ii,:] = -n_tot/rho*m*D[ii,:] + np.sum(m**2/rho**2*n_tot*n*D[ii,:])
    return D_Bruno

# <codecell>
#def q_mp_simple(m, p, n, Q, B=0, kT=0):
def q_mp_simple(*args):
    ''' Parameters q^mp as in Devoto, Phys. Fluids, 10, 2105 (1967)
        input may be: (m, p, n, Q) or (m, p, n, Q, B, kT)
        n: must contain species number density, with:
            n[0]:electron density
            n[1]:H+ density
            n[2]:H density
        analogous for m...
    '''
    if (len(args)!=4 and len(args)!=6):
        raise ValueError("Wrong number of input values!")
    m = args[0]; p = args[1]; n = args[2]; Q = args[3]
    
    # Since q^mp are symmetric (q^mp = q^pm) in the theory, I only define
    # the q^mp values where m >= p
    if (m>p):
        raise ValueError("[q_mp_simple] it must be p>=m, but m={},p={}".format(m,p))
    
    sumAQ = n[0]**2*np.sqrt(2) * np.sum(A[:,:,m,p]*Q[:,:,0,0])
    # [Opt] I could do this in a smarter way
    sumBQ = 0
    for j in range(1,len(n)):
        sumBQ += n[j]*np.sum(np.sum(B[0,:,m,p]*Q[0,:,0,j]))
    sumBQ *= n[0]

    q_mp = sumAQ + sumBQ
    
#    # I add effect of magnetic field (see Devoto,PhysFluids 11,448(1968), eq. 6)
#    if (len(args)==6):
#        mag_field = args[5]
#        kT = args[6]
#        omega_el = gau.qe*mag_field/m[0]
#        q_mp += -1j*n[0]*omega_el*np.sqrt(2*np.pi*m[0]/kT)*2/np.sqrt(np.pi)*(p+3/2)fattoriale!!!/pfattoriale * delta m,p
    if len(args)==6:
        raise ValueError("Magentic field correction not implemented yet!")
        
    return q_mp

def q_mp_complete(m, p, n_i, m_i, Q):
    ''' Parameters q^mp as in Devoto, Phys. Fluids, 9, 12305 (1966)
        n_i: must contain species number density
        m_i: must contain species atomic mass(only monatomic scpecies are allowed!)
    '''
    
    if not(len(n_i)==len(m_i)==Q.shape[2]==Q.shape[3]):
        raise ValueError("length of n_i must be same as length of m_i, and Q.shape[2] and Q.shape[3]")
    
    Nspecs = len(n_i)
    idx = np.arange(Nspecs)
    q = np.zeros((Nspecs,Nspecs))
    
    # I define conveniently a difference of Kronecker deltas:
    # ij_m_jl = delta_ij - delta_jl
    # l MUST BE A NUMPY ARRAY
    ij_m_jl = lambda i,j,l: ((i==j)-(j==l).astype(int))
    
    # I define conveniently a sum of Kronecker deltas:
    # ij_m_jl = delta_ij + delta_jl
    # l MUST BE A NUMPY ARRAY
    ij_p_jl = lambda i,j,l: ((i==j)+(j==l).astype(int))
    
    # q^00_ij
    if m==p==0:
        sqm = np.sqrt(m_i)
        for ii,jj in it.product(idx,idx):
            q[ii,jj] = 8*np.sum(n_i * sqm[ii]/np.sqrt(m_i[ii]+m_i)
                               *Q[0,0,ii,idx] * (n_i[ii]*sqm/sqm[jj]*ij_m_jl(ii,jj,idx)
                                               -n_i[jj]*sqm*sqm[jj]/m_i[ii]*(1-(ii==idx).astype(int))))
    # q^01_ij and q^10_ij
    elif (m==0 and p==1) or (m==1 and p==0):
        m32 = m_i**(3/2)
        for ii,jj in it.product(idx,idx):
            q[ii,jj] = 8*n_i[ii]*(m32[ii]/m32[jj]) \
                       *np.sum(n_i*m32/(m_i[ii]+m_i)**(3/2)
                              *(5/2*Q[0,0,ii,idx] - 3*Q[0,1,ii,idx])*ij_m_jl(ii,jj,idx))
            if (m==1 and p==0):  # q^10_ij
                q[ii,jj] *= m_i[jj]/m_i[ii]
    # q^11_ij
    elif (m==1 and p==1):
        m2 = m_i**2
        for ii,jj in it.product(idx,idx):
            q[ii,jj] = 8*n_i[ii]*(m_i[ii]/m_i[jj])**(3/2) \
                       *np.sum(n_i * np.sqrt(m_i)/(m_i[ii]+m_i)**(5/2)
                              *(ij_m_jl(ii,jj,idx) * (5/4*(6*m2[jj] + 5*m2)*Q[0,0,ii,idx]
                                                     -15*m2*Q[0,1,ii,idx]
                                                     +12*m2*Q[0,2,ii,idx])
                               +ij_p_jl(ii,jj,idx) * 4*m_i[jj]*m_i*Q[1,1,ii,idx]))
    # q^02_ij and q^20_ij
    elif (m==0 and p==2) or (m==2 and p==0):
        m52 = m_i**(5/2)
        for ii,jj in it.product(idx,idx):
            q[ii,jj] = 8*n_i[ii]*m52[ii]/m52[jj] \
                       *np.sum(n_i*m52/(m_i[ii]+m_i)**(5/2)
                              *ij_m_jl(ii,jj,idx)*(35/8*Q[0,0,ii,idx]
                                                  -21/2*Q[0,1,ii,idx]
                                                  +6*Q[0,2,ii,idx]))
            if (m==2 and p==0):  # q^20_ij
                q[ii,jj] *= (m_i[jj]/m_i[ii])**2
    # q^12_ij and q^21_ij
    elif (m==1 and p==2) or (m==2 and p==1):
        m2 = m_i**2
        for ii,jj in it.product(idx,idx):
            q[ii,jj] = 8*n_i[ii]*(m_i[ii]/m_i[jj])**(5/2) \
                       *np.sum(n_i*m_i**(3/2)/(m_i[ii]+m_i)**(7/2)
                              *(ij_m_jl(ii,jj,idx)*(35/16*(12*m2[jj]+5*m2)*Q[0,0,ii,idx]
                                                   -63/2*(m2[jj]+5/4*m2)*Q[0,1,ii,idx]
                                                   +57*m2*Q[0,2,ii,idx]
                                                   -30*m2*Q[0,3,ii,idx])
                                +ij_p_jl(ii,jj,idx)*(14*m_i[jj]*m_i*Q[1,1,ii,idx]
                                                    -16*m_i[jj]*m_i*Q[1,2,ii,idx])))
            if (m==2 and p==1):  # q^21_ij
                q[ii,jj] *= m_i[jj]/m_i[ii]
    # q^22_ij
    elif (m==2 and p==2):
        m2 = m_i**2
        m3 = m_i**3
        m4 = m_i**4
        for ii,jj in it.product(idx,idx):
            q[ii,jj] = 8*n_i[ii]*(m_i[ii]/m_i[jj])**(5/2) \
                       *np.sum(n_i*np.sqrt(m_i)/(m_i[ii]+m_i)**(9/2)
                              *(ij_m_jl(ii,jj,idx)*(35/64*(40*m4[jj]+168*m2[jj]*m2 +35*m4)*Q[0,0,ii,idx]
                                                   -21/8*m2*(84*m2[jj]+35*m2)*Q[0,1,ii,idx]
                                                   +3/2*m2*(108*m2[jj]+133*m2)*Q[0,2,ii,idx]
                                                   -210*m4*Q[0,3,ii,idx]
                                                   +90*m4*Q[0,4,ii,idx]
                                                   +24*m2[jj]*m2*Q[2,2,ii,idx])
                                +ij_p_jl(ii,jj,idx)*(7*m_i[jj]*m_i*(4*m2[jj]+7*m2)*Q[1,1,ii,idx]
                                                    -112*m_i[jj]*m3*Q[1,2,ii,idx]
                                                    +80*m_i[jj]*m3*Q[1,3,ii,idx])))
    
    return q

def makeQ(kT, ne):
    ''' To fill the required collision integrals for a mixture of
        e-,H+,H, for the third approx of the thermal conductivity,
        and resisitivity.
    '''
    T = kT/gau.kB
    # since: max(l)=2, max(s)=5, max(i)=max(j)=3, then:
    CollInt = np.zeros((3,5,3,3))
    
    # 'which' is a tuple of tuples containg all the (l,s,i,j) places
    # that must be filled inside Q
#    # FOR USE WITH SIMPLIFIED TREATEMENT
#    which = (  # e-e
#             (2,2,1,1),(2,3,1,1),(2,4,1,1),
#               # e-H+
#             (1,1,1,2),(1,2,1,2),(1,3,1,2),(1,4,1,2),(1,5,1,2),
#               # e-H
#             (1,1,1,3),(1,2,1,3),(1,3,1,3),(1,4,1,3),(1,5,1,3),
#               # H+-H+
#             (1,1,2,2),(1,2,2,2),(1,3,2,2),(1,4,2,2),(1,5,2,2),
#             (2,2,2,2),(2,3,2,2),(2,4,2,2),
#             (3,3,2,2),
#               # H+-H
#             (1,1,2,3),(1,2,2,3),(1,3,2,3),(1,4,2,3),(1,5,2,3),
#             (2,2,2,3),(2,3,2,3),(2,4,2,3),
#             (3,3,2,3),
#               # H-H
#             (1,1,3,3),(1,2,3,3),(1,3,3,3),(1,4,3,3),(1,5,3,3),
#             (2,2,3,3),(2,3,3,3),(2,4,3,3),
#             (3,3,3,3))
    which = (  # e-e
             (1,1,1,1),(1,2,1,1),(1,3,1,1),(1,4,1,1),(1,5,1,1),
             (2,2,1,1),(2,3,1,1),(2,4,1,1),
             (3,3,1,1),
               # e-H+
             (1,1,1,2),(1,2,1,2),(1,3,1,2),(1,4,1,2),(1,5,1,2),
             (2,2,1,2),(2,3,1,2),(2,4,1,2),
             (3,3,1,2),
               # e-H
             (1,1,1,3),(1,2,1,3),(1,3,1,3),(1,4,1,3),(1,5,1,3),
             (2,2,1,3),(2,3,1,3),(2,4,1,3),
             (3,3,1,3),
               # H+-H+
             (1,1,2,2),(1,2,2,2),(1,3,2,2),(1,4,2,2),(1,5,2,2),
             (2,2,2,2),(2,3,2,2),(2,4,2,2),
             (3,3,2,2),
               # H+-H
             (1,1,2,3),(1,2,2,3),(1,3,2,3),(1,4,2,3),(1,5,2,3),
             (2,2,2,3),(2,3,2,3),(2,4,2,3),
             (3,3,2,3),
               # H-H
             (1,1,3,3),(1,2,3,3),(1,3,3,3),(1,4,3,3),(1,5,3,3),
             (2,2,3,3),(2,3,3,3),(2,4,3,3),
             (3,3,3,3))

    for x in which:
        CollInt[x[0]-1,x[1]-1,x[2]-1,x[3]-1] = Qhat(x[0],x[1],x[2],x[3],
                                                    T,
                                                    ne)
    # Now I simmetrize with respect to the species
    # It seems that the collision integral have the time reversal property,
    # thus I can exchange safely the species
    for ll,ss in it.product(range(CollInt.shape[0]),range(CollInt.shape[1])):
        for ii,jj in it.combinations(range(CollInt.shape[2]),2) :
            if (CollInt[ll,ss,ii,jj] != CollInt[ll,ss,jj,ii]):
                if CollInt[ll,ss,ii,jj]==0e0:
                    CollInt[ll,ss,ii,jj] = CollInt[ll,ss,jj,ii]
                elif CollInt[ll,ss,jj,ii]==0e0:
                    CollInt[ll,ss,jj,ii] = CollInt[ll,ss,ii,jj]
                else:
                    Warning("Collision integral CollInt[{},{},{},{}] is not simmetric!".format(ll,ss,ii,jj))
    
    return CollInt

def Qhat(l, s, i, j, T, ne):
    ''' Collision integral (= πσ²Ω^(l,s)*)
        i,j: 1=electron
             2=H+
             3=H
    '''
    if (i==1 and j==2):
        return Q_eHp(l, s, T, ne)
    elif (i==j==1):
        return Q_ee(l, s, T, ne)
    elif (i==j==2):
        return Q_HpHp(l, s, T, ne)
    elif (i==j==3):
        return Q_HH(l, s, T)
    elif (i==1 and j==3):
        return Q_eH(l, s, T)
    elif ( (i==2 and j==3) or (i==3 and j==2)):
        return Q_HpH(l, s, T)
    else:
        raise ValueError("Wrong choice for i,j")

def Q_eHp(l, s, T, ne):
    ''' Collision integral for e-H+ as in Hahn, Mason, Phys. Fluids 14, 278 (1971),
        (see formula 51 of that paper) (it's the same as used by Bruno, "Transport
        properties...Mars..", ESA STR-256 2010)
    '''
    
    kT = T*gau.kB
    # Electrons Debye length
    lDeb_el = np.sqrt(kT/(4*np.pi*ne*gau.qe**2))
    # Ions+Electrons Debye length
#    lDeb_el = np.sqrt(kT/(4*np.pi*2*ne*gau.qe**2))
    
    phi0 = gau.qe**2/lDeb_el
    
    Tstar = kT/phi0
    Nl = ( 1 - (1+(-1)**l)/(2*(l+1)) )**-1
    
    As = 0.0
    for ss in range(1,s):
        As += s**-1
    
    Cl = 0
    if (l%2):
        for ll in range(1, l+1, 2):
            Cl += ll**-1
        Cl -= 1/(2*l)
    else:
        for ll in range(1, l, 2):
            Cl += ll**-1
    
    # I hope this interpretation of the paper by Hahn is correct
    # (my doubt is: is gamma already the euler_gamma constant?
    #  or is ln(gamma) the euler_gamma constant, as stated in the paper?)
    gamma2 = np.exp(2*np.euler_gamma)
    
    OmegaStar = (Tstar**-2) * \
                (l*Nl/(s*(s+1))) * \
                np.log((4*Tstar/gamma2)*np.exp(As-Cl)+1)
    Qls = np.pi*lDeb_el**2*OmegaStar
    return Qls

def Q_ee(l, s, T, ne):
    return Q_eHp(l, s, T, ne)

def Q_HpHp(l, s, T, ne):
    return Q_eHp(l, s, T, ne)

def Q_eH(l, s, T):
    '''Collision integrals for interaction e-H as in Bruno, Phys. Plasmas
    17, 112315 (2010).
    See eq. 19 and table VI'''
    
    g = g_eH[(l,s)]
    x = np.log(T)
    x1 = (x-g[0])/g[1]
    
    Sigma2OmegaStar = g[2] * x**g[4] * np.exp(x1) / (np.exp(x1)+np.exp(-x1)) + \
                      g[5]*np.exp(-((x-g[6])/g[7])**2) + g[3]
    # Convert to gaussian units
#   Sigma2OmegaStar *= 1e4 # assuming it was in meters
    Sigma2OmegaStar *= 1e-16 # assuming it was in Angstrom
    
    return Sigma2OmegaStar*np.pi

def Q_HH(l, s, T):
    ''' Collision integrals for interaction H-H as in Bruno, Phys. Plasmas
    17, 112315 (2010).
    See eq. 11 and table I'''

    a = a_HH[(l,s)]
    x = np.log(T)
    x1 = (x-a[2])/a[3]
    x2 = (x-a[5])/a[6]
    
    Sigma2OmegaStar = (a[0]+a[1]*x)*np.exp(x1)/(np.exp(x1) + np.exp(-x1)) + \
                       a[4]*np.exp(x2)/(np.exp(x2) + np.exp(-x2))
    # Convert to gaussian units
#    Sigma2OmegaStar *= 1e4
    Sigma2OmegaStar *= 1e-16 # assuming it was in Angstrom
    
    return Sigma2OmegaStar*np.pi

def Q_HpH(l, s, T):
    ''' Collision integrals for interaction H+-H in Bruno, Phys. Plasmas
    17, 112315 (2010).
    See eq. 11 and table I'''
    
    a = a_HpH[(l,s)]
    x = np.log(T)
    x1 = (x-a[2])/a[3]
    x2 = (x-a[5])/a[6]
    
    Sigma2OmegaStar = (a[0]+a[1]*x)*np.exp(x1)/(np.exp(x1) + np.exp(-x1)) + \
                      a[4]*np.exp(x2)/(np.exp(x2) + np.exp(-x2))
    # Odd-l therms must include the effect of inelastic collisions.
    # See eq.16 of Bruno's paper
    # I do not include l=3, because also Bruno does not (he does not provide
    # the required parameters for l=3, at the moment I have no idea why)
    if l==1 :
        d = d_HHp[(l,s)]
        Sigma2OmegaStar_chex = d[0] + d[1]*x + d[2]*(x**2)
        Sigma2OmegaStar = np.sqrt(Sigma2OmegaStar**2 + Sigma2OmegaStar_chex**2)
    # Convert to gaussian units
#    Sigma2OmegaStar *= 1e4
    Sigma2OmegaStar *= 1e-16 # assuming it was in Angstrom

    return Sigma2OmegaStar*np.pi
    

# <codecell>
# For testing purposes
if __name__ == "__main__":
#    rho = 3.16116e-7
    kT = 20000*gau.kB
    rho = 3.16116e-7
#    rho = 3.28e-5
#    kT = 24000*gau.kB
    z = iz.ionizationSaha(rho, kT)
    prs = (1+z)*rho/(gau.mp+gau.me)*kT
    
    ne = iz.elec_dens(rho,z)
    n_i = np.array([ne, ne, rho/(gau.mp+gau.me)-ne])
    m_i = np.array([gau.me, gau.mp, gau.mp+gau.me])
    
    Q = makeQ(kT, ne)
#    q00_h = q_mp_complete(0, 0, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    q01_h = q_mp_complete(0, 1, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    q10_h = q_mp_complete(1, 0, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    q11_h = q_mp_complete(1, 1, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    q12_h = q_mp_complete(1, 2, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    q21_h = q_mp_complete(2, 1, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    q22_h = q_mp_complete(2, 2, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    q20_h = q_mp_complete(2, 0, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    q02_h = q_mp_complete(0, 2, n_i[1:], m_i[1:], Qh[:,:,1:,1:])
#    
    q00 = q_mp_simple(0, 0, n_i, Q)
    
    D_heavy = diffu(kT, rho, n_i[1:], m_i[1:], Q[:,:,1:,1:])
    DT_heavy = thermDiffu(kT, n_i[1:], m_i[1:], Q[:,:,1:,1:])
    
    D_all = diffu(kT, rho, n_i, m_i, Q)
    DT_all = thermDiffu(kT, n_i, m_i, Q)
    
    DA, DAT = diffuAmbANDthermDiffuAmb(kT, rho, n_i, m_i, Q)
    
    lprime = lambdaPrime(kT, n_i[1:], m_i[1:], Q[:,:,1:,1:])
    
    k_tot = thermCond_norm(z, kT, rho)
    
    eta2 = elRes_norm_2(z,rho,kT)
    eta3 = elRes_norm_3(z,rho,kT)
    
    print("-------------")
    print("rho={}, T={}, z={}".format(rho, kT/gau.kB, z))
    print("prs={} atm".format(prs/1e6))
    print("-------------")

    # With z=0.935347, rho=3.16116e-7, T=20000, correct eta is eta=9.6e-15,
    # according to Bruno,"...Jupiter...", Table VIII
    print("eta:{:g}".format(elRes_norm(z,rho,kT)))
    print("eta2={}".format(eta2))
    print("eta3={}".format(eta3))

    print("Q^(1,1)_eH={}".format(Q_eH(1, 1, 9000)))
    print("Q^(1,1)_eH={}".format(Q_eH(1, 1, 18000)))
    print("Q^(2,2)_HH={}".format(Q_HH(2, 2, 9000)))
    print("Q^(2,2)_HH={}".format(Q_HH(2, 2, 18000)))
    print("Q^(2,2)_ee={}".format(Q_ee(2, 2, 9000, 1e16)))
    
    # With z=0.935347, rho=3.16116e-7, T=20000, correct kappa is kappa=2.979e5,
    # according to Bruno,"...Jupiter...", Table VIII
    print("k:{:g}".format(thermCond_e_norm(z,rho,kT)))
    
#    print("q00={}".format(q00))
#    print("q00_h={}".format(q00_h))
#    print("q01_h={}".format(q01_h))
#    print("q10_h={}".format(q10_h))
#    print("q11_h={}".format(q11_h))
#    print("q21_h={}".format(q21_h))
#    print("q12_h={}".format(q12_h))
#    print("q22_h={}".format(q22_h))
#    print("q02_h={}".format(q02_h))
#    print("q20_h={}".format(q20_h))
    
    print("D_heavy={}".format(D_heavy))
    print("DT_heavy={}".format(DT_heavy))
    print("D_all={}".format(D_all))
    
    print("lambda'={}".format(lprime))
    print("k_tot={}".format(k_tot))
    
    print("D_all, Bruno's convention:\n{}".format(D_Devoto2Bruno(n_i,m_i,D_all)))
    
    # In Devoto, J Plas Phys 2,4,617-631(1968),Table 7 DA(1atm,20000K)=579cm²/s
    print("DA={}".format(DA))
    print("DATl={}".format(DAT))