# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:36:57 2017

@author: aeber
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 08 16:58:05 2017

@author: aeber
"""
#shows a plot comparing the probability of being inside the field of view
#for all the upcoming telescopes listed below
#also shows a histogram of stars
#both as a function of x_in (1e6)
import matplotlib.pyplot as plt
import pickle
import math
import Math
import csv

l = .76e-6
scope_L = "LUVOIR"
D_L = 11.7
IWA_L = 3.
OWA_L = 30.
scope_J = "JWST"
D_J = 6.5
IWA_J = 11.3
OWA_J = 840.
scope_W = "WFIRST"
D_W = 2.4
IWA_W = 3.
OWA_W = 10. 
scope_H = "HabEx"
D_H = 8.
IWA_H = 3.
OWA_H = 20. 
file_ = "hygdata_v3.csv"


if __name__ == "__main__":
    s, p_s = pickle.load(open("p_s_norm.p", "rb"))
    d = [] #distance in pc
    x = []
    X = []
    p_inside_L = []
    p_inside_J = []
    p_inside_W = []
    p_inside_H = []
    #stars x and d
    with open(file_, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        i = 0
        for row in reader:
            if i > 1:
                d_ = float(row[9])*Math.pc 
                root_lum = math.sqrt(float(row[33]))
                x_ = float(d_)/float(Math.au*root_lum*1.68)
                if d_ < 30.*Math.pc and x_ < .005e9:
                    d.append(d_/Math.pc)
                    x.append(x_/1e6)
            i += 1
    N = sum(p_s)
    x_sorted = []
    for i in range(len(x)):
        x_sorted.append(x[i])
    x_sorted.sort()
    for x_i in x_sorted:
        x_ = x_i*1e6
        p_L, p_J, p_W, p_H = 0, 0, 0, 0
        for i in range(len(s)):
            if s[i] > IWA_L*x_*l/float(D_L) and s[i] < OWA_L*x_*l/float(D_L):
                p_L += p_s[i]
            if s[i] > IWA_J*x_*l/float(D_J) and s[i] < OWA_J*x_*l/float(D_J):
                p_J += p_s[i]
            if s[i] > IWA_W*x_*l/float(D_W) and s[i] < OWA_W*x_*l/float(D_W):
                p_W += p_s[i]
            if s[i] > IWA_H*x_*l/float(D_H) and s[i] < OWA_H*x_*l/float(D_H):
                p_H += p_s[i]
        p_inside_L.append(p_L/N)
        p_inside_J.append(p_J/N)
        p_inside_W.append(p_W/N)
        p_inside_H.append(p_H/N)
    print "Luvoir: ", sum(p_inside_L)
    print "JWST: ", sum(p_inside_J)
    print "WFIRST: ", sum(p_inside_W)
    print "HabEx: ", sum(p_inside_H)
    
    fig_, ax1_ = plt.subplots()
    plt.title("Sensitivity by telescope")
    ax1_.hist(x, bins = 50, color = 'k')
    ax1_.set_xlabel("$x_{intrinsic}$ (1e6)")
    ax1_.set_ylabel("counts", color = "k")
    ax1_.tick_params('y', colors = 'k')
    ax2_ = ax1_.twinx()
    ax2_.plot(x_sorted, p_inside_L, 'r-', label = "LUVOIR")
    ax2_.plot(x_sorted, p_inside_J, 'b-', label = "JWST")
    ax2_.plot(x_sorted, p_inside_W, 'g-', label = "WFIRST")
    ax2_.plot(x_sorted, p_inside_H, 'm-', label = "HabEx")
    ax2_.set_ylabel("probability planet is inside angle", color = 'k')
    ax2_.tick_params('y', colors = 'k')
    plt.legend()
    fig_.tight_layout()
    plt.show()

    