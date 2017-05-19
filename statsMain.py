# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 14:10:01 2016

@author: aeber
"""

import stats
import time
import manipulateAll as mA
import matplotlib.pyplot as plt
import reflectionSignal1_1 as rS
direct = "C:/Users/aeber/Desktop/Python/NASA"
import sys
sys.path.insert(0, direct)
import Math

#returns a planet object with O content fraction x
def O2vary(x, o0):
    return mA.EarthObj(Ofrac = x, T = o0.T, S = o0.S, d = o0.d, R = o0.R)

def H2Ovary(x, o0):
    return mA.EarthObj(H2Ofrac = x, T = o0.T, S = o0.S, d = o0.d, R = o0.R)

def gasVary(x, y, o0):
    return mA.EarthObj(Ofrac = x, H2Ofrac = y, T = o0.T, S = o0.S, d = o0.d, R = o0.R)

def Rvary(x, o0):
    return mA.EarthObj(Rfrac = x, T = o0.T, S = o0.S, d = o0.d , R = o0.R)

def Tvary(x):
    return mA.EarthObj(T = x, S = 1e-5, d = 6.)

def dvary(x):
    return mA.EarthObj(T = 1e5, S = 1e-5, d = x)

def ResVary(x):
    return mA.EarthObj(T = 1e4, S = 1e-5, d = 6., R = x)

#returns the likelihood function for P(data|param[i]) where data is randomly
#generated for given wavelengths l (in nm), with width dl (m) and planet object o0
#searches over parameters in given params, paramfunc should take a parameter
#as an argument and return a planet object for that parameter
def searchParameterSpace(l, dl, params, paramfunc, o0, detector = rS.QE_EMCCD):
    L = []
    signal = mA.findSignal(o0, l, dl, detector = detector)[1]
    noise = mA.findNoise(o0, l, dl, S = o0.S, detector = detector)[1]
    data = stats.poissonPick(Math.arrayAdd(signal, noise))
    expected = []
    for i in range(len(params)):
        oi = paramfunc(params[i], o0)
        signal_ = mA.findSignal(oi, l, dl, detector = detector)[1]
        expected = Math.arrayAdd(signal_, noise)
        L.append(stats.lnlikelyRatio(data, expected))
    return L


#returns the likelihood surface for P(data|param[0][i], param[1][j]) where data
#is randomly generated for wavelength l (in nm) with width dl (m) for planet object o0
#searches over params in two dimensions params[0] and params[1]
#paramfunc should take as an argument two parameter values params[0][i] and params[1][j]
def searchParameterSpace2D(l, dl, params, paramfunc, o0, detector = rS.QE_HgCdTe):
    L, X, Y = [], [], []
    Lmin, i_min, j_min = False, 0, 0
    signal = mA.findSignal(o0, l, dl, detector = detector)[1]
    noise = mA.findNoise(o0, l, dl, S = o0.S, detector = detector)[1]
    data = stats.poissonPick(Math.arrayAdd(signal, noise))
    for i in range(len(params[0])):
        print i
        B0 = params[0][i]
        L_i = []
        X_i = []
        Y_i = []
        for j in range(len(params[1])):
            B1 = params[1][j]
            oij = paramfunc(B0, B1, o0)
            signal_ = mA.findSignal(oij, l, dl, detector = detector)[1]
            expected = Math.arrayAdd(signal_, noise)
            L_ = stats.lnlikelyRatio(data,expected)
            if L_<Lmin or not(Lmin):
                Lmin, i_min, j_min = L_, B0, B1
            L_i.append(L_)
            X_i.append(B0)
            Y_i.append(B1)
        L.append(L_i)
        X.append(X_i)
        Y.append(Y_i)
    return L, X, Y, Lmin, i_min, j_min


#returns the lower and upper bounds on 67% and 95% confidence intervals as a function
#of given observation parameters (paramsObs)
def searchObsParameterSpace(l, dl, paramsObs, paramfuncObs, params, paramfunc, detector = rS.QE_EMCCD):
    lower1 = []
    upper1 = []
    lower2 = []
    upper2 = []
    start = time.time()
    for i in range(len(paramsObs)):
        print "running MC for x = ", paramsObs[i]
        print "Time Remaining: ", Math.timeRemaining(time.time() - start, i, len(paramsObs))
        start = time.time()
        o0 = paramfuncObs(paramsObs[i])
        L = testEnsemble(l, dl, params, paramfunc, o0, n = 20, detector = rS.QE_EMCCD)
        left_1, right_1 = findClosest2Values(L, min(L) + 1)
        left_2, right_2 = findClosest2Values(L, min(L) + 2)
        lower1.append(params[left_1])
        upper1.append(params[right_1])
        lower2.append(params[left_2])
        upper2.append(params[right_2])
    return lower1, upper1, lower2, upper2


#returns the value of the liklihood ration function for the average of n
#searches of parameter space
def testEnsemble(l, dl, params, paramfunc, o0, n = 100, detector = rS.QE_EMCCD):
    L = []
    L_bar = []
    for i in range(n):
        print i
        L.append(searchParameterSpace(l, dl, params, paramfunc, o0, detector))
    for i in range(len(L[0])):
        L_i = 0
        for j in range(n):
            L_i += L[j][i]
        L_bar.append(float(L_i)/n)
    return L_bar


#tests an ensemble in two dimensions instead of just one
def testEnsemble2D(l, dl, params, paramfunc, o0, n = 10, detector = rS.QE_HgCdTe):
    L = []
    X = []
    Y = []
    min_ = 0
    i_ = 0
    j_ = 0
    L_bar = Math.matrix(len(params[0]), len(params[1]))
    start = time.time()
    for i in range(n):
        print "running MC sim:", i
        print "Time remaining: ", Math.timeRemaining(time.time() - start, i, n)
        start = time.time()
        results = searchParameterSpace2D(l, dl, params, paramfunc, o0, detector = detector)
        L.append(results[0])
        if results[3] < min_ or i == 0:
            min_ = results[3]
            i_ = results[4]
            j_ = results[5]
    for i in range(len(L[0])):
        for j in range(len(L)):
            for k in range(n):
                L_bar[j][i] += L[j][i][k]
            L_bar[j][i] /= float(n)
    X = results[1]
    Y = results[2]
    return L_bar, X, Y, min_, i_, j_


#returns the indices of the two values in matrix A that are closest to val
def findClosest2Values(A, val):
    left_val, left = A[0], 0
    right_val, right = A[1], 1
    for i in range(2, len(A)):
        dc = abs(A[i] - val)
        dr = abs(right_val - val)
        dl = abs(left_val - val)
        if dc < dr or dc < dl:
            if dr < dl:
                left_val, left = right_val, right
            right_val, right = A[i], i
    return left, right


if __name__ == "__main__":
    l = range(280, 4000, 10)
    dl = 1e-8
    x = []
    iterations = 200
    res = 50.
    for i in range(iterations):
        x.append((i + 1)/res)
    xObs = []
    for i in range(30,150,4):
        xObs.append(i)
    left1, right1, left2, right2 = searchObsParameterSpace(l, dl, xObs, ResVary, x, O2vary)
    plt.figure(1)
    plt.plot(xObs, left1, "r")
    plt.plot(xObs, right1, "r")
    plt.plot(xObs, left2, "b")
    plt.plot(xObs, right2, "b")
    plt.show()
    w = []
    for i in range(len(left1)):
        w.append(abs(left2[i] - right2[i]))
    plt.figure(2)
    plt.plot(xObs, w)
    plt.show()
    