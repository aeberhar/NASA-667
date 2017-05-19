# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 10:53:30 2017

@author: aeber
"""
#calculates and plots the probability distribution for a planet seperation from the star
#projected in the observational plane
#order is dr, de, da
import matplotlib.pyplot as plt
import math
import time
import Math

res = 75 #should be kept around 75 to avoid getting to close to function singularities

T_star = 1.0 #effective temperature of target star
a_in = 1.0/1.68
a_out = 1.0
T_2 = (T_star/Math.T_sun)**2
a_c = ((a_out + a_in)/2.)*T_2
a_max, a_min = a_out*T_2, a_in*T_2
rho = 2.*a_max*a_min/float(a_max + a_min)


if __name__ == "__main__":
    da, ds = float(a_max - a_min)/res, float(a_max)/res
    s, s_norm, a, p_s = [], [] , [], []
    for i in range(res - 1):
        p_s.append(0)
        s.append(ds*i + ds)
        s_norm.append(s[-1]/float(a_max))
    for i in range(res - 1):
        a.append(a_min + i*da + da)
    start_ = time.time()
    start = start_
    for y in range(len(a)): #loop for p(a)
        a_ = a[y]
        C, e = T_2/a_, []
        if a_ < ((a_out + a_in)/2.)*T_2:
            e_max = min(1. - a_in*C, .99)
        else:
            e_max = max(a_out*C - 1.,0.)
        de = e_max/res
        for i in range(res - 1):
            e.append(de*i + de)
        for e_ in e: #loop for p(e):
            r_min, dr, r = a_*(1. - e_), 2.*a_*e_/res, []
            for i in range(res - 1):
                r.append(r_min + dr*i + dr)
            for r_ in r: #loop for p(r|e,a)
                for x in range(len(s)): #loop for p(s|r)
                    s_ = s[x]
                    if s_ < r_:
                        dB = da*de*dr
                        num = math.sqrt(((r_*e_)**2 + 2.*r_*a_*(1. - e_**2) - r_**2))
                        denom = math.sqrt((((r_*e_)**2 - (a_*(1. - e_**2) - r_)**2)*(2./r_ - 1./a_)))
                        p_r_ae = float(num)/denom
                        p_s_r = (s_/float(r_))**2 /abs((a_**(3./2.))*math.sqrt(r_**2 - s_**2)*math.pi**2)
                        p_i = p_r_ae*p_s_r*dB
                        p_s[x] += p_i 
        dT = abs(start - time.time())
        print abs(start_ - time.time())/(60.*60.), " hours elapsed"
        print dT*abs(len(a) - y - 1)/(60.*60.), " hours remaining\n"
        start = time.time()
    plt.figure(2)
    plt.ylabel("probability density")
    plt.xlabel("apparent separation")
    plt.plot(s_norm, p_s)
    plt.show()
    
    