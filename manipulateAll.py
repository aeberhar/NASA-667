# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 17:19:33 2016

@author: aeber
"""
#an auxiliary file that uses all the files so far to produce results

import blackbody as bb
import reflectionSignal1_1 as rS
import candidate1_1 as candidate
import matplotlib.pyplot as plt
direct = "C:/Users/aeber/Desktop/Python/NASA"
import sys
sys.path.insert(0, direct)
import Math
import Albedo
import numpy.random as np


#finds and returns the noise at the point of the planet due to the PSD binned
#for given planet object and wavelength of light matrix l in nm with intervals dl in m
#bins photons according to given detector with given starlight supression factor (<1)
def findNoise(planet, l, dl, detector = rS.QE_EMCCD, S = 1):
    a, r, d, D, R = planet.a, planet.r_p, planet.d, planet.D, planet.R
    data_star = bb.readData(planet.L_file)
    L = []
    for i in range(len(l)):
        l_ = l[i]*1e-9
        L_star = [bb.data2I(l_, dl, data_star[0], data_star[1])[0], 0]
        L.append(rS.diffraction(a, r, l_, L_star, d = d, D = D)[0])
    F = rS.Lum2Photons(l, L, D, d)[0]
    Math.scaleMatrix(F, S*planet.T)
    return rS.photons2Bins(l, F, R, detector)    


#finds and returns the stellar luminosity of a star at given wavelengths l (in nm)
#for given planet object and wavelength intervals dl (in m)
#units of luminosity are W/m^2
def findStellarSpectrum(planet, l , dl, watts = False):
    data_ = False
    if planet.L_file:
        data_ = True
        data_star = bb.readData(planet.L_file)
    L = []
    for i in range(len(l)):
        l_ = l[i]*1e-9
        if data_:
            L.append(bb.data2I(l_, dl, data_star[0], data_star[1])[0])
        else:
            L.append(bb.I_(l_, planet.T_eff, dl)[0])
    if watts:
        Math.scaleMatrix(L, 4*Math.pi*planet.r_star[0]**2)
    return L


#find the ideal habital zone of a planet
#this is defined at the position where the 
def getHZ(o):
    L = bb.findL(o.T_eff, o.r_star)
    L_ratio = Math.math.sqrt(L[0]/Math.L_solar)
    return L_ratio*Math.au
        

#finds and returns the signal from the planet in binned photon fluxes
#for given planet object and wavelength of light matrix l with intervals dl
#bins photons according to given detector
def findSignal(planet, l, dl, detector = rS.QE_EMCCD):
    L = []
    data_star = bb.readData(planet.L_file)
    d, D, R = planet.d, planet.D, planet.R
    for i in range(len(l)):
        l_ = l[i]*1e-9
        A = Albedo.findA(l_,planet, .8, .3)[0]
        L_star = [bb.data2I(l_, dl, data_star[0], data_star[1])[0], 0]
        L.append(rS.reflect(A, [300,0],l_, dl, L_star, planet, total = 0)[0])
    F = rS.Lum2Photons(l, L, D, d)[0]
    Math.scaleMatrix(F, planet.T)
    return rS.photons2Bins(l, F, R, detector)   
        

#creates a generic earth like planet object with given parameters
def EarthObj(Ofrac = 1.0, H2Ofrac = 1.0, D = [10*Math.pc, 0], T = 100, d = 10, S = 1, Rfrac = 1.0, R = 100):
    o = candidate.ToI(True)
    o.N_str[1][1] *= Ofrac
    o.N_str[2][1] *= Ofrac
    o.N_tr[1][1] *= Ofrac
    o.N_tr[2][1] *= H2Ofrac
    o.r_p[0] *= Rfrac
    o.r_p[1] *= Rfrac
    o.R = R
    o.D = D
    o.T = T
    o.d = d
    o.S = S
    return o


#creates an exoplanet object with given parameter and randomly generates other variables ir rand is true
def exoPlanet(o0, Ofrac = 1.0, H2Ofrac = 1.0, rand = True):
    o = candidate.ToI(True)
    o.d = o0.d
    o.S = o0.S
    o.R = o0.R
    o.T = o0.T
    o.L_file = o0.L_file
    o.N_str[1][1] *= Ofrac
    o.N_str[2][1] *= Ofrac
    o.N_tr[1][1] *= Ofrac
    o.N_tr[2][1] *= H2Ofrac
    if o0.r_p[1] > 0 and rand:
        r_p = np.normal(o0.r_p[0], o0.r_p[1])
        while(r_p < 0):
            r_p = np.normal(o0.r_p[0], o0.r_p[1])
        o.r_p[0] = r_p
    else:
        o.r_p[0] = o0.r_p[0]
    if o0.a[1] > 0 and rand:
        a = np.normal(o0.a[0], o0.a[1])
        while(a < 0):
            a = np.normal(o0.a[0], o0.a[1])
        o.a[0] = a
    else:
        o.a[0] = o0.a[0]
    if o0.m[1] > 0 and rand:
        m = np.normal(o0.m[0], o0.m[1])
        while(m < 0):
            m = np.normal(o0.m[0], o0.m[1])
        o.m[0] = m
    else:
        o.m[0] = o0.m[0]
    if o0.T_eff[1] > 0 and rand:
        T_eff = np.normal(o0.T_eff[0], o0.T_eff[1])
        while(T_eff < 0):
            T_eff = np.normal(o0.T_eff[0], o0.T_eff[1])
        o.T_eff[0] = T_eff
    else:
        o.T_eff[0] = o0.T_eff[0]
    if o0.r_star[1] > 0 and rand:
        r_star = np.normal(o0.r_star[0], o0.r_star[1])
        while(r_star < 0):
            r_star = np.normal(o0.r_star[0], o0.r_star[1])
        o.r_star[0] = r_star
    else:
        o.r_star[0] = o0.r_star[0]
    if o0.D[1] > 0 and rand:
        D = np.normal(o0.D[0], o0.D[1])
        while(D < 0):
            D = np.normal(o0.D[0], o0.D[1])
        o.D[0] = D
    else:
        o.D[0] = o0.D[0]
    return o