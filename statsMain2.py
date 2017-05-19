# -*- codi, ng: utf-8 -*-
"""
Created on Wed Aug 03 09:42:35 2016

@author: aeber
"""
#this is an updated versio of statsMain, instead of Earth-like planets around
#sun-like stars I am modeling the nearest stars. Gaussian picks uncertain 
#parameters

import stats
import pickle
import time
import manipulateAll as mA
import candidate1_1 as candidate
import matplotlib.pyplot as plt
import reflectionSignal1_1 as rS
direct = "C:/Users/aeber/Desktop/newComputer/Python/NASA/667/"
import sys
sys.path.insert(0, direct)
import Math


def loadPlanet(i, name):
    inFile = open(name, "rb")
    planets = pickle.load(inFile)
    return planets[i]


#returns a planet that is Earth-like in the HZ of a star with given effective temperature (K)
#stellar radius and Kurusc model file, at distance D (m)
def modelPlanet(T_eff, r_star, L_file, D):
    o = candidate.ToI(True)
    o.T_eff = T_eff
    o.r_star = r_star
    o.L_file = L_file
    o.a = [mA.getHZ(o),.1*mA.getHZ(o)]
    o.L_file = L_file
    o.D = D
    o.S = 1e-5
    o.T = 1e3
    o.d = 6.
    o.R = 2.
    print "Time:", o.T
    print "d:", o.D
    print "D:", o.d
    print "S:", o.S
    print "Temp:", T_eff
    print "r:", r_star
    return o


#returns a planet object with O content fraction x
def O2vary(x, o0):
    return mA.exoPlanet(o0, Ofrac = x, rand = False)

def H2Ovary(x, o0):
    return mA.exoPlanet(o0, H2Ofrac = x, rand = False)


#returns the likelihood function for P(data|param[i]) where data is randomly
#generated for given wavelengths l (in nm), with width dl (m) and planet object o0
#searches over parameters in given params, paramfunc should take a parameter
#as an argument and return a planet object for that parameter
def searchParameterSpace(l, dl, params, paramfunc, o0, detector = rS.QE_EMCCD, binary = False):
    L = []
    o_rand = mA.exoPlanet(o0, rand = True)
    signal = mA.findSignal(o_rand, l, dl, detector = detector)[1]
    noise = mA.findNoise(o_rand, l, dl, S = o0.S, detector = detector)[1]
    if binary:
        noise1 = mA.findNoise(binary, l, dl, S = o0.S, detector = detector)[1]
        signal1 = mA.findSignal(binary, l, dl, detector = detector)[1]
        noise = Math.arrayAdd(noise, noise1)
        signal = Math.arrayAdd(signal, signal1)
    data = stats.poissonPick(Math.arrayAdd(signal, noise))
    expected = []
    for i in range(len(params)):
        oi = paramfunc(params[i], o0)
        signal_ = mA.findSignal(oi, l, dl, detector = detector)[1]
        if binary:
            oi_ = paramfunc(params[i], binary)
            signal_ = Math.arrayAdd(mA.findSignal(oi_, l, dl, detector = detector)[1], signal_)
        expected = Math.arrayAdd(signal_, noise)
        L.append(stats.lnlikelyRatio(data, expected))
    return L


#returns the value of the liklihood ration function for the average of n
#searches of parameter space
def testEnsemble(l, dl, params, paramfunc, o0, n = 100, detector = rS.QE_EMCCD, binary = False):
    L = []
    L_bar = []
    for i in range(n):
        print i
        L.append(searchParameterSpace(l, dl, params, paramfunc, o0, detector, binary))
    for i in range(len(L[0])):
        L_i = 0
        for j in range(n):
            L_i += L[j][i]
        L_bar.append(float(L_i)/n)
    return L_bar


def findClosest(A, val):
    i = 0
    d = abs(val - A[0])
    for i_ in range(len(A)):
        d_ = abs(val - A[i_])
        if d_ < d:
            d = d_
            i = i_
    return i


#prints out the confidence interval bounds
def findCI(L, x, val):
    mid = L.index(min(L))
    left = L[0:mid]
    right = L[mid:]
    i_left = findClosest(left, val)
    i_right = findClosest(right,val)
    return x[i_left], x[mid + i_right]
    

if __name__ == "__main__":
    #the "real" planet to be modeled
    o0 = modelPlanet([3042,117], [.141*Math.r_solar, 0.007*Math.r_solar], "kp00_3500.ascii", [4.246*Math.ly, .006*Math.ly])
    l = range(280, 4000, 10)
    dl = 1e-8
    x = []
    iterations = 200
    res = 50.
    for i in range(iterations):
        x.append((i + 1)/res)
    L = testEnsemble(l, dl, x, H2Ovary, o0, n = 1000, detector = rS.JWST_NIRCAM)
    L_min = min(L)
    plt.figure(1)
    plt.plot(x, L)
    plt.plot([x[0], x[-1]], [L_min + 1, L_min + 1])
    plt.plot([x[0], x[-1]], [L_min + 3.84, L_min + 3.84])
    plt.show()
    print "95% CI:", findCI(L, x, L_min + 3.84)
    print "Time:", o0.T
    print "d:", o0.D
    print "D:", o0.d
    print "S:", o0.S
    print "Temp:", o0.T_eff
    print "r:", o0.r_star
    