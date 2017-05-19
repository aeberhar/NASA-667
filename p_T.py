# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:06:27 2017

@author: aeber
"""
#calculates the transmission pdf for a given star and wavelength
import Math
import pickle
import matplotlib.pyplot as plt
import random
import math
import time
import numpy as np

trials = 999
d = 12.*Math.ly #distance to star,m
D = 11.7 #apeture diameter, m
l = .76e-6 #wavelength of light being observed in, m
OWA = 30*float(l*d)/D
IWA = 3.3*float(l*d)/D
L = .52*Math.L_sun
L_ = math.sqrt(L/float(Math.L_sun))
a_max = 1.68*Math.au*L_

#heagonal shear
def T_calc(x_, y_):
    return (math.sin(math.pi*D*y_/float(9.*l))*math.sin(math.pi*D*float(math.sqrt(3)*x_ + y_)/float(18.*l)))**2

#parallel shear
def T_calc_2(x_, y_):
    return (math.sin(math.pi*D*y_/float(9.*l)))**4
    

if __name__ == "__main__":
    s, p_s = pickle.load(open("p_s_norm.p", "rb"))
    T = []
    start_ = time.time()
    start = start_
    zeroes = 0
    T_xy, X, Y = [], [], []
    T_xy2 = []
    dX = 40.*l/float(D*5000)
    dx_ = 40./float(5000.)
    for i in range(5000):
        x_ = -20.*l/float(D) + dX*i
        _x_ = -20. + dx_*i
        X.append(_x_)
        Y.append(_x_)
        temp, temp2 = [], []
        for j in range(5000):
            y_ = -20.*l/float(D) + dX*j
            if math.sqrt(x_**2 + y_**2) > OWA/d:
                temp.append(0)
                temp2.append(0)
            else:
                temp.append(T_calc(x_, y_))
                temp2.append(T_calc_2(x_, y_))
        T_xy.append(temp)
        T_xy2.append(temp2)
    for i in range(trials):
        phi_ = random.uniform(0, 2*math.pi)
        s_ = Math.pick(s, p_s)*Math.au
        if s_ > OWA or s_ < IWA:
            T.append(0)
            zeroes += 1
        else:
            x_, y_ = s_*math.cos(phi_), s_*math.sin(phi_)
            T_ = T_calc(x_,y_)
            if T_**2 == 0:
                print "wow"
                zeroes += 1
            T.append(T_**2)
        if i%(trials/10) == 0:
            dT = abs(start - time.time())
            print int(((i+1)/float(trials))*100), " percent complete"
            print abs(start_ - time.time())/(60.*60.), " hours elapsed"
            print dT*abs(trials - i - 1)/(60.*60.*(trials/10)), " hours remaining\n"
            start = time.time()
    print np.average(T), "+/-", 1.96*np.std(T)/math.sqrt(trials)
    z_ = zeroes/float(trials)
    print "chance of 0 transmission: " , z_, "+/-", math.sqrt((1-z_)*z_)/math.sqrt(trials) 
    plt.figure(1)
    plt.xlabel("Transmission")
    plt.ylabel("counts")
    a = plt.hist(T, bins = 1000)
    plt.show()
    plt.figure(2)
    plt.xlabel("seperation ($\lambda$/D)")
    plt.ylabel("seperation ($\lambda$/D)")
    plt.contourf(X, Y, T_xy)
    plt.colorbar()
    plt.show()
    plt.figure(3)
    plt.xlabel("seperation ($\lambda$/D)")
    plt.ylabel("seperation ($\lambda$/D)")
    plt.contourf(X, Y, T_xy2)
    plt.colorbar()
    plt.show()