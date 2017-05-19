# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 11:20:29 2017

@author: aeber
"""
#looks at the transmission probability distribution function for a specific star
import Math
import pickle
import matplotlib.pyplot as plt
import random
import math
from astropy.io import fits
import time
import numpy as np

trials = 999999
d = 12.*Math.ly #distance to star,m, this star is Tau Ceti
D = 11.7 #apeture diameter, m
L = .52*Math.L_sun #luminosity of target star
L_ = math.sqrt(L/float(Math.L_sun)) 
a_max = 1.68*Math.au*L_ #maximum semi major axis



if __name__ == "__main__":
    s, p_s = pickle.load(open("p_s_norm.p", "rb"))
    start_ = time.time()
    start = start_
    l, dl = Math.fill_array(.1e-6, 5.0e-6, 200) #wavelength of light being observed in, m
    ex_T, z_T = [], [] #expectation value for T, chance of getting 0
    i = 0
    for l_ in l: #for each wavelength
        T = []
        zeroes = 0
        OWA = 20*float(l_*d)/D
        IWA = 4.0*float(l_*d)/D
        for y in range(trials): #run trials
            phi_ = random.uniform(0, 2*math.pi)
            s_ = Math.pick(s, p_s)*a_max
            if s_ > OWA or s_ < IWA:
                T.append(0)
                zeroes += 1
            else:
                x_, y_ = s_*math.cos(phi_), s_*math.sin(phi_)
                T_ = math.sin(math.pi*D*y_/float(9.*l_))*math.sin(math.pi*D*float(math.sqrt(3)*x_ + y_)/float(18.*l_))
                if T_**2 == 0:
                    print "wow"
                    zeroes += 1
                T.append(T_**2)
        ex_T.append(np.average(T))
        z_T.append(zeroes/float(trials))
        if i%(len(l)/10) == 0:
            dT = abs(start - time.time())
            print int(((i+1)/float(len(l)))*100), " percent complete"
            print abs(start_ - time.time())/(60.*60.), " hours elapsed"
            print dT*abs(len(l) - i - 1)/(60.*60.*(len(l)/10)), " hours remaining\n"
            start = time.time()
        i += 1
    l_microns = []
    for i in range(len(l)):
        l_microns.append(l[i]*1e6)
    plt.figure(1)
    plt.xlabel("wavelength (microns)")
    plt.ylabel("average Transmission")
    plt.plot(l_microns, ex_T)
    plt.show()
    plt.figure(2)
    plt.xlabel("wavelength (microns)")
    plt.ylabel("probability of 0 transmission")
    plt.plot(l_microns, z_T)
    plt.show()