# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:44:24 2017

@author: aeber
"""
#plots the integration space in the e-a and a-r planes for the p(s) function
import matplotlib.pyplot as plt
import Math
from astropy.io import fits
file_ = "4RING_2SHEAR_union.fits"

a_max = 1.0
a_min = 1.0/1.68
res = 500

da = float(a_max - a_min)/res
a_c = float(a_max + a_min)/2.
rho = 2.*a_max*a_min/float(a_max + a_min)
e_max = float(a_max - a_min)/(a_max + a_min)


if __name__ == "__main__":
    a, e = [], []
    A = fits.open(file_)
    ap = A[0].data
    for i in range(res + 1):
        a_ = da*i + a_min
        a.append(a_)
        if a_ <= a_c:
            e.append(1. - a_min/a_)
        else:
            e.append(a_max/a_ - 1.)
    plt.figure(1)
    plt.xlabel("semi major axis (fraction of max)")
    plt.ylabel("three space seperation (fraction of max)")
    plt.plot([0, 2], [max(e), max(e)], color = 'b')
    plt.plot([a_min, a_min], [0, max(e) + 1], color = 'k')
    plt.plot([a_max, a_max], [0, max(e) + 1], color = 'k')
    plt.plot(a, e, color = 'r')
    plt.show()
    plt.figure(2)
    plt.xlabel("semi major axis (fraction of max)")
    plt.ylabel("eccentricity")
    plt.plot([a_min, a_min], [a_min, a_max], color = 'k')
    plt.plot([a_max, a_max], [a_min, a_max], color = 'k')
    plt.plot([a_min, a_max], [a_min, a_min], color = 'k')
    plt.plot([a_min, a_max], [a_max, a_max], color = 'k')
    plt.plot([a_min, a_c], [rho, a_max], color = 'r')
    plt.plot([a_c, a_max], [a_min, rho])
    plt.show()
    plt.figure(3)
    plt.contourf(ap)
    plt.show()