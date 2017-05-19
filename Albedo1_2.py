# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:49:34 2016

@author: aeber
"""

#finds the abledo of a planet, entirely model based, results are exact
#model approximates atmosphere as a series of two dimensional planes
#stratsphere and troposphere layers act as semi-opaque rayleigh scatterers
#cloud and surface act as semi-mirrors
#auxillary file

import Math
import gasData as gD

#finds and returns the Rayleigh scattering optical depth (dimensionless), exact
def findR(N, t, l):
    return ((337.55)*(1.12e-29)**2)*(N[1]*t/(l**4))


#finds and returns the cross section for this gas at this particular wavelength
#if the cross section data does not include this particlar wavelength we return the closest
def getSig(l, gas):
    cross, l_ = gD.cross[gas][1], gD.cross[gas][0]
    i = Math.binarySearchClose(l, l_)
    return cross[i]


#returns reflection,transmission, attenuation coefficients for a specific atmospheric layer, exact
#given atmospheric composition matrix N, with N[gas][gas name, atoms/m^3], exact
#layer thickness in meters, exact; wavelength of light l in m, exact
def findFrT(N, t, l):
    tau = 0
    R = 0
    for i in range(len(N)):
        sig = getSig(l, N[i][0]) #cross section for this particular gas
        tau += sig*N[i][1]*t #opticl thickness of this layer due to absorption
        R += findR(N[i], t, l) #rayleigh scattering optical thickness
    abs_ = Math.e**(-1.0*tau) #fraction not absorbed
    r = (1.0 - Math.e**(-1.0*R)) #light reflected by layer due to Rayleigh scattering
    return abs_*r, abs_*(1-r), (1 - abs_) #reflect, transmitted to surface, absorbed fractions
    

#returns cloud albedo at a given wavelength, exact
#using approximation detailed in Kokhanovsky, 2002 "Optical properties of terrestial clouds"
def getCloudA(l):
    sigma_o = 1.0 #this is sigma_0 * 10 m
    dsig = 1282
    l_ = l**(2.0/3.0)
    sigma = sigma_o*(1.0 + l_*dsig)
    return (1.0 - Math.e**(-1.0 * sigma))


#finds and returns reflected and trasmitted coefficients,exact
#given flux at top of atmosphere and uncertainties for a given wavelength
#l is given wavelength is exact (given in meters not nm!!!)
#planet object
#cloud albedo, cloud fraction, surface albedo (A_c, f_c, A_s) are exact
def findA(l, planet, f_c, A_s):
    try:
        N_str, t_str, N_tr, t_tr = planet.N_str, planet.t_str, planet.N_trop, planet.t_trop
    except AttributeError:
        N_str, t_str, N_tr, t_tr = planet.N_str, planet.t_str, planet.N_tr, planet.t_tr
    #stratosphere
    r_str, T, abs_str = findFrT(N_str, t_str, l)
    #clouds
    A_c = getCloudA(l)
    r_c = T*f_c*abs_str*A_c #light reflected by clouds that escapes atmosphere
    T *= (1. - r_c) #fraction that is transmitted through clouds
    #troposphere
    args_tr = findFrT(N_tr, t_tr, l)
    r_tr = args_tr[0]*T*(1 - abs_str) #reflected fraction, reflection coef*transmitted to that layer*not aborbed by stratosphere
    T *= args_tr[1] #fraction that is transmitted through troposphere, reaches sea level
    #surface
    r_s = T*A_s*(1. - args_tr[2])*(1. - abs_str) #transmitted*albedo*absorbed by trop and strat
    r = r_s + r_tr + r_c + r_str
    T_ = 1. - r #fraction that goes into heating planet
    return r, T_, T