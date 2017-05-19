# -*- coding: utf-8 -*-
"""
Created on Tue Aug 09 10:52:55 2016

@author: aeber
"""
#this code is supposed to make confidence belts and it doesn't but im not sure why
#all my parameters always converge to 1 instead of whatever they're actual value is

import stats
import pickle
import time
import manipulateAll as mA
import candidate1_1 as candidate
import matplotlib.pyplot as plt
import reflectionSignal1_1 as rS
direct = "C:/Users/aeber/Desktop/Python/NASA"
import sys
sys.path.insert(0, direct)
import Math


def loadPlanet(i, name):
    inFile = open(name, "rb")
    planets = pickle.load(inFile)
    return planets[i]


#returns a planet that is Earth-like in the HZ of a star with given effective temperature (K)
#stellar radius and Kurusc model file, at distance D (m)
def modelPlanet(T_eff, r_star, L_file, D, Ofrac = 1.0):
    o = candidate.ToI(True)
    o.N_str[1][1] *= Ofrac
    o.N_str[2][1] *= Ofrac
    o.N_tr[1][1] *= Ofrac
    o.T_eff = T_eff
    o.r_star = r_star
    o.L_file = L_file
    o.a = [mA.getHZ(o),.1*mA.getHZ(o)]
    o.L_file = L_file
    o.D = D
    o.S = 1e-5
    o.T = 1e3
    o.d = 6.
    o.R = 100.
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
def searchParameterSpace(l, dl, params, paramfunc, o0, detector = rS.QE_EMCCD, binary = False, Ofrac = 1.0, H2Ofrac = 1.0):
    L = []
    o_rand = mA.exoPlanet(o0, Ofrac = Ofrac, H2Ofrac = H2Ofrac, rand = True)
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
            signal_ = Math.arrayAdd(mA.findSIgnal(oi_, l, dl, detector = detector)[1], signal_)
        expected = Math.arrayAdd(signal_, noise)
        L.append(stats.lnlikelyRatio(data, expected))
    return L


#returns the value of the liklihood ration function for the average of n
#searches of parameter space
def testEnsemble(l, dl, params, paramfunc, o0, n = 100, detector = rS.QE_EMCCD, binary = False, Ofrac = 1.0, H2Ofrac = 1.0):
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
    l = range(280, 4000, 10)
    dl = 1e-8
    y = []
    iterations_ = 16
    res_ = 8.
    for i in range(iterations_):
        y.append((i + 1)/res_)
    x = []
    iterations = 160
    res = 40.
    for i in range(iterations):
        x.append((i + 1)/res)
    left = []
    right = []
    start = time.time()
    for i in range(len(y)):
        print "running MC sim:", i
        print "Time remaining: ", Math.timeRemaining(time.time() - start, i, len(y))
        start = time.time()
        o0 = modelPlanet([3970,0], [.605*Math.r_solar, 0.002*Math.r_solar], "kp00_4000.ascii", [3.587*Math.pc, .01*Math.pc], Ofrac = y[i])
        L = testEnsemble(l, dl, x, O2vary, o0, Ofrac = y)
        L_min = min(L)
        bounds = findCI(L, x, L_min + 3.84)
        left.append(bounds[0])
        right.append(bounds[1])
    plt.figure(1)
    for i in range(len(left)):
        plt.plot([left[i], right[i]], [y[i], y[i]])
    plt.show()
    