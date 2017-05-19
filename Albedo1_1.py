# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 14:46:48 2016

@author: aeber
"""

import Math
import gasData as gD
import blackbody as bb


#finds and returns the average mass per atom of your atmospheric layer, exact
def getM(N):
    M = 0.
    norm = 0.
    for i in range(len(N)):
        gas = N[i][0]
        factor = N[i][1]
        norm += factor
        M += gD.M[gas]*factor
    return M/norm


#reursive algorithm
#returns the reflection coefficient of a layer
#given height above the bottom of that layer z, exact; gas mixing matrix N
#differential bin width dz, exact; scaling height H with uncertainty
#number density at bottom of layer N0, exact; recursive depth i
def reflection(z, N, dz, H, N0, i, l):
    if i == 0:
        return r_(z, dz, H, N0, l), a_(z, N, dz, H, N0, l)
    else:
        r_m, a_m = reflection(z + dz, N, dz, H, N0, i - 1, l)
        r_n = r_(z, dz, H, N0, l)
        a_n = a_(z, N, dz, H, N0, l)
        return doubleSemiMirrorR(r_m, r_n, a_m, a_n), doubleSemiMirrorAbs(r_m, r_n, a_m, a_n)


#returns the reflective coefficient for light rayleigh scattered by a layer
#given height of section z, exact; differential bin width, exact
#scaling height H with uncertainty; density at bottom of layer N0, exact
def r_(z, dz, H, N0, l):
    arg = getRSig(l)*-.5*dz
    n = N_(z, H, N0)
    R = Math.scale(n[0], n[1], arg)
    b = Math.exp(R[0], R[1])
    return 1. - b[0], b[1]


#returns the absortion coefficient for light absorbed by a layer
#given height of section z, exact; differential bin width, exact
#scaling height H with uncertainty; density at bottom of layer N0, exact
def a_(z, N, dz, H, N0, l):
    arg = getSig(N, l)*-1*dz
    n = N_(z, H, N0)
    tau = Math.scale(n[0], n[1], arg)
    b = Math.exp(tau[0], tau[1])
    return 1. - b[0], b[1]


#finds and returns the cross section for this gas at this particular wavelength
#if the cross section data does not include this particlar wavelength we return the closest
def getGasSig(l, gas):
    cross, l_ = gD.cross[gas][1], gD.cross[gas][0]
    i = Math.binarySearchClose(l, l_)
    return cross[i]

#returns the integrated cross section given a gas matrix that contains mixing ratios
def getSig(N, l):
    sig = 0.
    norm = 0.
    for i in range(len(N)):
        norm += N[i][1]
        sig += getGasSig(l, N[i][0])*N[i][1]
    return sig/norm


#calculates rayleigh scattering cross section
#follows work of Bodhaine, B., Wood, N., Dutton, E., Slusser, J. 1999, Journal
#of Atmospheric and Oceanic Technology, 16, 1854.
def getRSig(l):
    return ((337.55)*(1.12e-29)**2)*(1.0/(l**4))
    

#returns the density at a given height above sea level (z) and scale height (H)
#z is exact, H and returned density should have uncertainties
#N0 is the number density of air at sea level
def N_(z, H, N0):
    arg = Math.divide(z, 0, H[0], H[1])
    exp = Math.exp(arg[0], arg[1])
    return Math.scale(exp[0], exp[1], N0)


#returns the effective reflective coefficient for two parrallel mirrors
#given reflective coefficients of first and second mirrors (r1 and r2) with uncertainty
#and absorption coefficient for each mirror (a1 and a2) with uncertainty
def doubleSemiMirrorR(r1, r2, a1, a2):
    T = Math.add(1. - r1[0], r1[1], -1*a1[0], a1[1])
    T2 = Math.square(T[0], T[1])
    num = Math.multiply(r2[0], r2[1], T2[0], T2[1])
    denom = Math.multiply(r1[0], r1[1], r2[0], r2[1])
    inf = Math.divide(num[0], num[1], 1. - denom[0], denom[1])
    return Math.add(r1[0], r1[1], inf[0], inf[1])


#returns the effective absorption coefficient for two parrallel mirrrors
#given reflective coefficients of first and second mirrors (r1 and r2) with uncertainty
#and absorption coefficient for each mirror (a1 and a2) with uncertainty
def doubleSemiMirrorAbs(r1, r2, a1, a2):
    T = Math.add(1. - r1[0], r1[1], -1*a1[0], a1[1])
    ar = Math.multiply(r2[0], r2[1], a1[0], a1[1])
    T_ = Math.add(ar[0], ar[1], a2[0], a2[1])
    num = Math.multiply(T_[0], T_[1], T[0], T[1])
    denom = Math.multiply(r1[0], r1[1], r2[0], r2[1])
    inf = Math.divide(num[0], num[1], 1. - denom[0], denom[1])
    return Math.add(a1[0], a1[1], inf[0], inf[1])


#given planet object (planet), stellar luminosity (L) with uncertainty
#average atomic mass of atoms in air in kg/atom (M), exact
#returns the scale height with uncertainty
def findH(planet, L, M):
    r, m, a = planet.r_p, planet.m, planet.a
    C = Math.k/(2.0*Math.G*(Math.steph_boltz*Math.pi)**(.25))/M
    a2 = Math.square(a[0], a[1])
    r2 = Math.square(r[0], r[1])
    F = Math.divide(L[0], L[1], a2[0], a2[1])
    F_4 = Math.sqrt(F[0], F[1], .25)
    Fmoment = Math.multiply(r2[0], r2[1], F_4[0], F_4[1])
    prop = Math.divide(Fmoment[0], Fmoment[1], m[0], m[1])
    return Math.scale(prop[0], prop[1], C)


#returns cloud albedo at a given wavelength, exact
#using approximation detailed in Kokhanovsky, 2002 "Optical properties of terrestial clouds"
def getCloudA(l):
    sigma_o = 1.0 #this is sigma_0 * 10 m
    dsig = 1282
    l_ = l**(2.0/3.0)
    sigma = sigma_o*(1.0 + l_*dsig)
    return (1.0 - Math.e**(-1.0 * sigma))


def findA(l, planet, f_c, A_s, N0, N_str, N_trop, i = 100):
    M = getM(N_trop) #average mass per atom of our gas
    L = bb.findL(planet.T, planet.r_star)
    H = findH(planet, L, M) #scale height in meters
    N0_str = N0/(Math.e**3)
    t_trop = 3.*H
    t_str = 4.*H
    dz_str = float(t_str)/i #differential width for stratosphere
    r_str, a_str = reflection(0, N_str, dz_str, H, N0_str, i, l) #effective reflection ofstratosphere
    r_c = [getCloudA(l)*f_c,0]
    r1 = doubleSemiMirrorR(r_str, r_c, a_str, [0,0]) #strat and clouds
    a1 = doubleSemiMirrorAbs(r_str, r_c, a_str, [0,0])
    dz_trop = float(t_trop)/i
    r_trop, a_trop = reflection(0, N_trop, dz_trop, H, N0, i, l)
    r2 = doubleSemiMirrorR(r1, r_trop, a1, a_trop)
    a2 = doubleSemiMirrorAbs(r1, r_trop, a1, a_trop)
    return doubleSemiMirrorR(r2, A_s, a2, [0,0])