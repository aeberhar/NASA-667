# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 16:59:52 2016

@author: aeber
"""

import candidate1_1 as candidate
import blackbody as bb
import reflectionSignal1_1 as rS
import matplotlib.pyplot as plt
direct = "C:/Users/aeber/Desktop/Python/NASA"
import sys
sys.path.insert(0, direct)
import Math
import Albedo
import stats

#sim stuff
l = range(280, 4000, 1) #wavelengths of light sim will include, (nm)
dl = 1.0e-9

def EarthObj(a = [Math.au, .1*Math.au], T_eff = [5770, 50], L_file = "kp00_5750.ascii", atm = True):
    o = candidate.ToI()
    o.a = a
    o.D = [10.0*Math.pc, Math.pc*.15]
    o.m = [Math.m_earth, .1*Math.m_earth]
    o.r_p = [Math.r_earth, .1*Math.r_earth]
    o.r_star = [Math.r_solar, .1*Math.r_solar]
    o.T_eff = T_eff
    o.L_file = L_file
    o.d = 10
    if atm:
        o.N_str = [["N2", 2.0e22], ["O3", 1.5e20], ["O2", 4.0e21]] #atmospheric gases matrix [gas][gas name, atoms/m^3], stratosphere, exact
        o.t_str = 30.0e3 #thickness meters, stratosphere, exact
        o.N_tr = [["N2", 2.0e25], ["O2", 5.0e24], ["H2O", 1.0e22]] #atmospheric gases matrix for troposphere, exact
        o.t_tr = 20.0e3
    return o


if __name__ == "__main__":
    planet0 = EarthObj()
    planet1 = EarthObj([.98*Math.au, 0],[6030,0],"kp00_6000.ascii")
    planet2 = EarthObj([1.63*Math.au, 0],[6030,0],"kp00_6000.ascii")
    planet3 = EarthObj(atm = False)
    data_solar = bb.readData(planet0.L_file)
    data_star = bb.readData(planet1.L_file)
    L0 = []
    L1 = []
    L2 = []
    L3 = []
    L_sun = []
    A_ = []
    A_noatm_ = []
    contrast = []
    for x in l:
        l_ = x*1e-9
        L_solar = [bb.data2I(l_, dl, data_solar[0], data_solar[1])[0], 0]
        L_sun.append(L_solar[0])
        L_star = [bb.data2I(l_, dl, data_star[0], data_star[1])[0], 0]
        A = Albedo.findA(l_,planet0, .8, .3)[0]
        A_noatm = Albedo.findA(l_, planet3, 0, .3)[0]
        A_.append(A)
        A_noatm_.append(A_noatm)
        L0_stuff =rS.reflect(A, [300,0],l_, dl, L_solar, planet0, total = 0, contrast = True) 
        contrast.append(L0_stuff[2][0])
        L0.append(L0_stuff[0][0])
        L1.append(rS.reflect(A, [300,0],l_, dl, L_star, planet1, total = 0)[0])
        L2.append(rS.reflect(A, [273,0],l_, dl, L_star, planet2, total = 0)[0])
        L3.append(rS.reflect(A_noatm, [300,0],l_,dl, L_solar, planet3, total = 0)[0])
    #luminosity of planets
    plt.figure(1)
    plt.plot(l,L0)
    plt.plot(l, L1)
    plt.plot(l ,L2)
    plt.plot(l ,L3)
    plt.show()
    #converted to photons
    plt.figure(2)
    contrast = rS.Lum2Photons(l, contrast, planet0.D, planet0.d)[0]
    F0 = rS.Lum2Photons(l, L0, planet0.D, planet0.d)[0]
    F1 = rS.Lum2Photons(l, L1, planet0.D, planet0.d)[0]
    F2 = rS.Lum2Photons(l, L2, planet0.D, planet0.d)[0]
    F3 = rS.Lum2Photons(l, L3, planet0.D, planet0.d)[0]
    plt.plot(l, F0)
    plt.plot(l, F1)
    plt.plot(l, F2)
    plt.plot(l, F3)
    plt.show()
    #convert photons to binned photons
    plt.figure(3)
    contrast = rS.photons2Bins(l, contrast, 100., rS.QE_EMCCD)
    N0 = rS.photons2Bins(l, F0, 100., rS.QE_EMCCD)
    N1 = rS.photons2Bins(l, F1, 100., rS.QE_EMCCD)
    N2 = rS.photons2Bins(l, F2, 100., rS.QE_EMCCD)
    N3 = rS.photons2Bins(l, F3, 100., rS.QE_EMCCD)
    plt.bar(N0[0], N0[1], color = 'red')
    plt.show()
    plt.figure(4)
    plt.bar(contrast[0], contrast[1])
    plt.show()