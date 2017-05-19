# -*- coding: utf-8 -*-
"""
Created on Wed Feb 08 16:58:05 2017

@author: aeber
"""
#makes plots showing the probability of being inside the field of view as a function of x
#and the position of the IWA and OWA compared to the on-sky pdf as a function of x
import matplotlib.pyplot as plt
import pickle
import math
import Math
import csv

l = .76e-6 #wavelength being observed at
scope = "LUVOIR" #telescope
D = 11.7 #aperture diameter
IWA = 3. #inner working angle
OWA = 30. #outer working angle
file_ = "hygdata_v3.csv" #csv file with information about stars within 30 pc


if __name__ == "__main__":
    s, p_s = pickle.load(open("p_s_norm.p", "rb"))
    inner = []
    outer = []
    d = []
    x = []
    X = []
    p_inside = []
    #inner and out working angle lines
    for i in range(len(s)):
        inner.append(s[i]/IWA)
        outer.append(s[i]/OWA)
    #stars x and d
    with open(file_, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        i = 0
        for row in reader:
            if i > 1:
                d_ = float(row[9])*Math.pc 
                root_lum = math.sqrt(float(row[33]))
                x_ = float(l*d_)/float(D*Math.au*root_lum*1.68)
                if d_ < 30.*Math.pc and x_ < 1./IWA:
                    d.append(d_/Math.pc)
                    x.append(x_)
            i += 1
    N = sum(p_s)
    x_sorted = []
    for i in range(len(x)):
        x_sorted.append(x[i])
    x_sorted.sort()
    for x_ in x_sorted:
        p_ = 0
        for i in range(len(s)):
            if s[i] > IWA*x_ and s[i] < OWA*x_:
                p_ += p_s[i]
        p_inside.append(p_/N)
    print "inner product of probability of being in field of view with star histogram: ", sum(p_inside)
        
    fig, ax1 = plt.subplots()
    ax1.plot(s,p_s, 'b-')
    ax1.set_xlabel("seperation (fraction of maximum)")
    ax1.set_ylabel('probability density', color = 'b')
    ax1.tick_params('y', colors = 'b')
    ax2 = ax1.twinx()
    ax2.plot(s, inner, 'r-')
    ax2.plot(s, outer, 'r-')
    ax2.set_ylabel('x', color = 'r')
    ax2.tick_params('y', colors = 'r')
    fig.tight_layout()
    plt.show()
    
    fig_, ax1_ = plt.subplots()
    plt.title(scope)
    ax1_.hist(x, bins = 100)
    ax1_.set_xlabel("x")
    ax1_.set_ylabel("counts", color = "b")
    ax1_.tick_params('y', colors = 'b')
    ax2_ = ax1_.twinx()
    ax2_.plot(x_sorted, p_inside, 'r-')
    ax2_.set_ylabel("probability planet is inside angle", color = 'r')
    ax2_.tick_params('y', colors = 'r')
    fig_.tight_layout()
    plt.show()

    