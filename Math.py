# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:39:59 2017

@author: aeber
"""
#a file containing physical constants and some useful functions

import random

#all constants are in SI (mks) units
T_sun = 5578. #effective temperature of sun, K
au = 1.496e11 #astronomical unit, m
pc = 3.086e16 #parsec, m
ly = 9.461e15 #light year, m
G = 6.67408e-11 #Newtons constant, m^3/(kg*s^2)
r_e = 6.371e6 #earth radius, m
r_sun = 695.7e6 #radius of sun, m
h = 6.626e-34 #planks constant, J*s
c = 3.0e8 #speed of light, m/s
k = 1.38064852e-23 #Boltzmann constant m^2kg /s^2K
L_sun = 3.828e26 #luminosity of sun, W
M_sun = 1.989e30 #mass of sun, kg
steph_boltz = 5.6703e-8 #stephan boltzmann constant, W/m^2K^4

#picks a random x given pdf y(x)
def pick(x, y):
    looking = True
    high = max(y)
    while(looking): 
        x_ = random.randint(0, len(x) - 1) #pick a random x value
        y_ = random.uniform(0, high) #pick a random y value
        if y_ < y[x_]: #if this point is inside the pdf then return that x value
            looking = False
            return x[x_]

#makes an array from min_ to max_ with res elements
def fill_array(min_, max_, res):
    A = []
    dx = float(max_ - min_)/res
    for i in range(res - 1):
        A.append(dx*i + dx + min_)
    return A, dx
    
    
#returns the index of val in sorted array A
#if val is not in A it returns the index of the val closest to val
def binarySearchClose(val, A):
    return binarySearchClosePartner(val, A, 0, len(A) - 1)

def binarySearchClosePartner(val, A, left, right): #"private partner"
    i_ = (right - left)/2 + left
    if A[i_] == val: #value is found
        return i_
    elif left >= right or i_ == 0 or i_ == len(A) - 1: #value is not in arrray
        if (i_ == 0 and A[i_] > val) or (i_ == len(A) - 1 and A[i_] < val): #end cases
            return i_
        elif val < A[i_]: #value is higher than current selection
            return i_ - int(abs(A[i_] - val) > abs(A[i_- 1] - val))
        elif val > A[i_]: #value is lower than current selection
            return i_ + int(abs(A[i_] - val) > abs(A[i_+ 1] - val))
    elif A[i_] < val: #selection is higher than value, search continues
        return binarySearchClosePartner(val, A, i_ + 1, right)
    else: #selection is lower than value, search continues
        return binarySearchClosePartner(val, A, left, i_ - 1)
        