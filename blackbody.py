# -*- coding: utf-8 -*-
"""
Created on Thu Jun 09 18:51:02 2016

@author: Andillio
"""
#this function can be used to make a blackbody spectrum
import astropy.io.ascii as ascii
direct = "C:/Users/aeber/Desktop/newComputer/Python/NASA/667"
import sys
sys.path.insert(0, direct)
import Math

#power per area for a given frequency (not angular), exact
#at a given temperature T, with uncetainty
#returns W/m^2 (i.e. steradians are integrated out)
def I(f, T, df):
    C = Math.pi*df*Math.h*(f**3)*2./float(Math.c**2)
    _T, s_T = Math.divide(1., 0, T[0], T[1])
    BE, sBE = Math.scale(_T, s_T, Math.h*f/Math.k)
    boltzF, s_ = Math.pow(Math.math.e, 0, BE, sBE)
    return C/(boltzF - 1.), C*s_

#same thing as other equation but for a given wavelength
#T given in K with uncertainties, dl is the integration bin width
#returns W/m^2 (i.e. steradians are integrated out)
def I_(l, T, dl):
    try:
        C = Math.pi*dl*Math.h*2.*(Math.c**2)/(l**5)
        _T, s_T = Math.divide(1., 0, T[0], T[1])
        BE, sBE = Math.scale(_T, s_T, (Math.h*Math.c)/(Math.k*l))
        boltzF, s_ = Math.pow(Math.math.e, 0, BE, sBE)
        return C/(boltzF - 1.), C*s_
    except ZeroDivisionError:
        print l, T, dl


#returns a luminosity of a star at given wavelength l and wavelength width dl
#finds closest value in array of available wavelengths wl, finds closest
#luminosity in L (units are W/m)
def data2I(l, dl, wl, L):
    i_ = Math.binarySearchClose(l, wl)
    return L[i_]*dl
    

#returns total luminosity of a star given its effective temperature T;
#given effective stellar radius r, with uncertainties
def findL(T, r):
    T4 = Math.square(T[0], T[1], 4)
    r2 = Math.square(r[0], r[1])
    A = Math.multiply(T4[0], T4[1], r2[0], r2[1])
    C = 4.*Math.pi*Math.steph_boltz
    return A[0]*C, A[1]*C


#returns a wavelength and luminosity per wavelength matrices
#given ascii data file name, constant c to convert wavelengths to m
#given constant C to convert luminosity to mks units
def readData(fileName, c = 1e-10, C = (1e7)):
    l, I = [], []
    data = ascii.read(fileName)
    for i in range(len(data)):
        l.append(data[i][0]*c)
        I.append(data[i][1]*C)
    return l, I
    