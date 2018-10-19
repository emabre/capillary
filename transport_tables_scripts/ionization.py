# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 17:37:37 2017

@author: ema
"""

import numpy as np
import constantsGAU_ema as gau
import scipy.optimize as opt
import importlib

importlib.reload(gau)

# Ionization energy for the single H atom in erg
chi = 13.597*gau.eV2erg
chi_diss = 4.476*gau.eV2erg
# chi = 15.398*gau.eV2erg

# Definition of free energy F corrected from in PR E 65 016407 (in the neutral gas free energy component
# you have to write 4V instead of V (reason: the hydrogen atom 1s ground state is four-fold degenerate
# (both the proton and electron have spin up or spin down, so its partition funciton should be
# multiplied by a factor 4 [Gilmore, "Saha's Equation and the Partition Function", pag 3])))
def F(z, chi_0, T_e, T_i, V, A, Epsilon_0, Nu_0):
    '''Function defined by eq 18 in PR E 65 016407, but corrected
    (in the neutral gas free energy component you have to write 4V instead of V)'''
    return z*chi_0 - z*T_e*( 1+np.log( (gau.me*T_e/(2*np.pi*gau.hbar**2))**(3/2) * (2*V/z) ) ) - z*T_i*( 1+np.log( (A*gau.mp*T_i/(2*np.pi*gau.hbar**2))**(3/2) * (2*V/z) ) ) - (1-z)*T_i*(1+np.log( (A*gau.mp*T_i/(2*np.pi*gau.hbar**2))**(3/2) * (4*V/(1-z)) )) + Epsilon_0/6 * (2+(Nu_0/V)**3 - 3*Nu_0/V )

# First derivative of F
def DF(z, chi_0, T_e, T_i, V):
    '''Derivative of F with respect to z (corrected from PR E 65 016407)'''
    return chi_0 - T_e*np.log((gau.me*T_e/(2*np.pi*gau.hbar**2))**(3/2) *2*V/z) + T_i*np.log(2*z/(1-z))

#
# BE CAREFUL!! TO THE UNITS OF THE DATA!
#
def ionizationBob(chi_0_loc, T_e_loc, T_i_loc, V_loc, tolerance):
    '''Computes ionization degree for a plasma (2 temeprature model, in principle it could be good also for species different from hydrogen, but I tested it only for hydrogen) on a certain point whose characteristics are given by the input parameters. Temperatures in erg, all the rest in usual cgs-gaussian units'''
    DF_z2energy = lambda z_loc: DF(z_loc, chi_0_loc, T_e_loc, T_i_loc, V_loc)
    z_min = 0+tolerance
    z_max = 1-tolerance
    if DF_z2energy(z_min)>=0 and DF_z2energy(z_max)>=0:
        z_loc = tolerance
    elif (DF_z2energy(z_min)<=0 and DF_z2energy(z_max)<=0):
        z_loc = 1-tolerance
    else:
        z_loc = opt.bisect(DF_z2energy, z_min, z_max, xtol=tolerance, maxiter=100)
    return z_loc

def ionizationSaha(rho, kT):
    '''Computes ionization degree for
    hydorgen according to Saha's
    model. Input: rho in g/cm^3,
    kT in erg'''
    n0 = rho/gau.mp  # total number of H atoms
    f = (2*np.pi*gau.me*kT)**(1.5)/(gau.h**3 * n0)*np.exp(-chi/kT)
    z = 0.5*(-f + np.sqrt(f**2+4*f))
    # print("chi:{}".format(chi))
    # print("f:{}".format(f))
    return z

def ionizDissSaha(rho, kT, ret_n0_H2=False):
    ''' Computes rigorously the ionization and dissociation degree for hydrogen,
        with Saha model.
        Returns a tuple: (dissociation, ionization)
        Definition of:
            (i)  Ionization (ys): ys is such that: ne = 2*ys*n0_H2,
                 where n0_H2 is the number density of Hydrogen H2 molecules which
                 may form at most (i.e.: twice the number density of hydrogen atoms
                 present, in whatever form, H+, H, H- ... whatever)
            (ii) Dissociation (xs): xs is such that: n_H2 = (1-xs)*n0_H2,
                 where n_H2 is the actual number density of H2 molecules present
                 at the specified conditions.
        Note that the two definisions imply that n_H = 2*(xs-ys)*n0_H2,
        where n_H is the number density of H atoms (non ionized) present at the
        specified conditions.
    '''
    
    # For convenience, I change the sign to the dissociation and ionization energy
    E1 = -chi_diss
    E2 = -chi
    
    n0_H = rho/(gau.me+gau.mp)  # total number of hydrogen atoms per cm³
    n0_H2 = 0.5*n0_H
    m_H = gau.mp + gau.me

    beta = 1/(kT)
    A = np.exp(beta*E1 - 3/2*np.log(-beta*E1) \
        + 3/2*np.log(2*np.pi*(0.5*m_H)*(-E1)/(n0_H2**(2/3)*gau.h**2)) )
    B = np.exp(beta*E2 - 3/2*np.log(-beta*E2) \
        + 3/2*np.log(2*np.pi*gau.me*(-E2)/(n0_H2**(2/3)*gau.h**2)) )
    
    # Ionization and Dissociation degree
    py = [16/(A*B), 0, 2, B, -B]
    ys = take_only_acceptable(np.roots(py))
    
    # xs = ys + 2*ys**2/(B*(1+1e-20))
    # xs = ys*(1 + 2*ys/B)
    # This expression for xs works better than the other two (for numerical reasons)
    xs = (-(A-8*ys) + np.sqrt((A-8*ys)**2 - 4*4*(4*ys**2 - A)) )/8
    
    if ret_n0_H2:
        return xs, ys, n0_H2
    else:
        return xs, ys

def take_only_acceptable(y):
    '''Returns only acceptable value for ionization degree among some values (input type: np.array).
    If none is found, it throws exception
    If more than one are acceptable but below a certain tolerance it returns the first one'''
    toll = 1e-4
    # I take only real numbers
    y_real = []
    for vv in y:
        y_maybe_real = np.real_if_close(vv)
        if not(np.iscomplexobj(y_maybe_real)):
            y_real.append(y_maybe_real)
    # No I select, among the real ones, the ones inside [0,1]
    if len(y_real)==0:
        raise ValueError("[find_only_acceptable]No real values:{}".format(y))
    y_ok = []
    for vv in y_real:
        try:
            if (vv>=0 and vv<=1):
                y_ok.append(vv)
        except:
            pass
    if len(y_ok)>1:
        if abs(min(y_ok)-max(y_ok))<toll:
            return y_ok[0]
        else:
            raise ValueError("[find_only_acceptable]Too many acceptable vaues:{}".format(y_ok))
    else:
        return float(y_ok[0])

def moleFract_Hfulldiss(z):
    '''Computes molar fractions of H, H+, e- for fully dissociated,
       partially ionized, hydrogen.
       Returns a np array, y:y[0]:H's molar fraction;
                             y[1]:H+'s molar fraction
                             y[2]:e-'s molar fraction'''
    return np.array([(1-z)/(1+z), z/(1+z), z/(1+z)])

def elec_dens(rho,z):
    return z * rho / (gau.mp+gau.me)
