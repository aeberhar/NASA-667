# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 20:09:49 2017

@author: aeber
"""
#runs a simulation in which each star is visited once and observed until a planet is found (or excluded) by each telescope
#outputs the constraints on eta that could be acheived given such a campaign
import matplotlib.pyplot as plt
import numpy as np
import math

eta = [.0001, .001, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.2, 1.5, 1.9, 2., 3.]
#xi values come from running makePlots or makePlots2
xi = [710.49, 53.31, 85.04] #number of stars we expect to see if we look at all and they all have 1
c = ['r', 'b', 'g'] #used for confidence interval colors
labels = ["LUVOIR", "JWST", "WFIRST"] #various telescopes
trials = 9999

left, right = [[],[],[]], [[],[],[]]
diff = [[],[],[]]

if __name__ == "__main__":
    for j in range(len(xi)): #telescope loop
        xi_ = xi[j]
        for eta_ in eta: #eta loop
            results = [] #results of estimations of true value of eta
            I_results = []
            for i in range(trials): 
                N_ = eta_*xi_
                N_i = np.random.poisson(N_) #number of planets found
                results.append(N_i/xi_)
                I_results.append(3.36/math.sqrt(N_))
            results.sort()
            """left[j].append(results[int(trials/4.)])
            right[j].append(results[int(3.*trials/4.)])
            diff[j].append(right[j][-1] - left[j][-1])"""
            mu = np.mean(results)
            s = np.std(results)
            left[j].append(mu - s*1.68)
            right[j].append(mu + s*1.68)
            #diff[j].append(np.mean(I_results))
            diff[j].append(3.36*s)
        print xi_, " done"
    for j in range(len(xi)): #plot each telescope results
        plt.figure(1)
        plt.ylabel("bounds on eta")
        plt.xlabel("true value of eta")
        plt.plot(eta, left[j], color = c[j], label = labels[j])
        plt.plot(eta, right[j], color = c[j])
        plt.legend()
        plt.show()
        plt.figure(2)
        plt.ylabel("width of confidence interval")
        plt.xlabel("true value of eta")
        plt.plot(eta, diff[j], color = c[j], label = labels[j])
        plt.legend()
        plt.show()