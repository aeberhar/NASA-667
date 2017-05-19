# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 11:47:49 2016

@author: aeber
"""
#contains functions that determine the reflection signal of a planet
#auxiliary file

import Math
import blackbody as bb

#calculates and returns reflected luminosity of a planet
#A is planet albedo at a given wavelength and is exact
#L is luminosity of host star a given wavelength and uncertainty
#a is orbital radius of planet and uncertainty
#r is planet radius and uncertainty
#face on fraction f, star light supression factor S (<1)
#total is boolean toggle for including starlight diffraction or not
#contrast is a boolean toggle for returning an extra value for the luminosity
#contrast at the point of the planet
def reflect(A, T_eff, l, dl, L_, planet, S = 1, f = .5, total = True, contrast = False):
    a, r, D, d = planet.a, planet.r_p, planet.D, planet.d
    r2 = Math.square(r[0], r[1])
    a2 = Math.square(a[0], a[1])
    r2 = Math.scale(r2[0], r2[1], Math.pi)
    a2 = Math.scale(a2[0], a2[1], Math.pi * 4.)
    L = Math.multiply(L_[0], L_[1], r2[0], r2[1])
    L = Math.scale(L[0], L[1], A*f)
    L = Math.divide(L[0], L[1], a2[0], a2[1])
    L_heat = findLHeat(T_eff, Math.scale(r2[0], r2[1], 4.0*f), l, dl)
    L_diff = diffraction(a, r, l, L_, d = d, D = D)
    L_diff = Math.scale(L_diff[0], L_diff[1], S)
    L_less = Math.add(L[0], L[1], L_heat[0], L_heat[1])
    L = L_less
    if total:
        L = Math.add(L[0], L[1], L_diff[0], L_diff[1])
    if contrast:
        return L, float(L_less[0])/L_diff[0], L_diff
    else:
        return L
    

#returns luminosity at a given wavelength l, temperature T
#and surface area a and integration bin width dl
def findLHeat(T, a, l, dl):
    F = bb.I_(l, T, dl)
    L = Math.multiply(F[0], F[1], a[0], a[1])
    return L


#returns the luminosity due to the star's diffraction through the optical
#equipment, given an orbital radius a, apeture diameter d (exact), 
#wavelength l (exact), distance to system D and luminosity of star L
#r is the planetary radius
def diffraction(a, r, l, L_, d = 10., D = [1.5*Math.pc, 0]):
    SMD_ = Math.scale(D[0], D[1], l/d) #smallest resolvable distance
    theta = Math.divide(a[0], a[1], SMD_[0], SMD_[1])
    P_ = PSF(theta)
    return Math.multiply(P_[0], P_[1], L_[0], L_[1])
    

#point spread function for stellar luminosity, (unitless)
#given theta value (defined as theta_seperation/theta_)
def PSF(theta):
    C = [0,0]
    C[0] = 1.0/(1.0 + ((Math.pi**4)/8.0)*theta[0]**3)
    C[1] = (3.0/8.0)*theta[1]*(Math.pi**4)*(theta[0]**2)*C[0]**2
    P = Math.scale(C[0], C[1], (Math.pi**(2./3.))*(Math.math.sqrt(27.)/16.))
    return P


#converts an array of luminosities (L, in W) at given wavelengths (l, in m)
#at a given distance with uncertainties (D), to a photon flux at the point
#of the detector (photons/s) and an array of uncerainties
def Lum2Photons(l, L, D, d = 10):
    F = []
    sF = []
    A = d**2
    D2 = Math.square(D[0], D[1])
    c = float(A)/(Math.h*Math.c*Math.pi*4.)
    C = Math.divide(c, 0, D2[0], D2[1])
    for i in range(len(l)):
        l_ = l[i]
        if l_ > 1:
            l_ *= 1e-9
        L_ = L[i]
        F_, s_ = Math.scale(C[0], C[1], l_*L_)
        F.append(F_)
        sF.append(s_)
    return F, sF


#bins photon flux in given array of photon fluxes (F) at wavelengths (l in nm), exact
#bins according to spectral resolution R, exact
#returns the left edge of each bin and the photon fluxes in each bin
def photons2Bins(l, F, R, e = lambda x:1.):
    Bins, lefts = [0], [l[0]]
    start, dl, j = l[0], float(l[0])/R, 0
    for i in range(len(l)):
        l_ = l[i]
        if l_ < (start + dl):
            Bins[j] += F[i]*e(l_*1e-9)
        else:
            start, dl, j = l_, float(l_)/R, j + 1
            Bins.append(F[i]*e(l_*1e-9))
            lefts.append(l_)
    return lefts, Bins


#returns quantum efficiency for EMCCD type detector 
def QE_EMCCD(l):
    C = .9*(2.*Math.math.sqrt(2.*Math.pi))*1e-7
    return C*Math.gauss(l, 2.0e-7, 6.0e-7)


def JWST_NIRCAM(l):
    if l > .6e-6 and l < 5e-6:
        return 1.0
    else:
        return 0.0


#returns the quantum efficiency for a HgCdTe detector
def QE_HgCdTe(l):
    return .7*(1. - Math.U(l - 1.1e-6)) + .8*Math.U(l - 1.1e-6)