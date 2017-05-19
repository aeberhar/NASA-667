# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 15:09:49 2017

@author: aeber
"""

#calculates how long it would take to exclude 95% of HZExs around a given star
#saves the list of star objects for a given telescope (if last line is not commented out)

import csv
import pickle
import Math
import math
import rand_orbits as ro
import matplotlib.pyplot as plt
import numpy as np
import time

#telescope stats
l, dl = .76e-6, .76e-7 #observed wavelength, and bandwidth
D = 8. #aperture diameter, m
scope = "HabEx" #the name of the mission, must be JWST, HabEx or LUVOIR
T = .1 #transmission
T_tot = 60.*60.*24.*365. #total mission lifetime, seconds
C = 1.0e-9 #contrast
IWA, OWA = 3., 20. #inner and outer working angles in l/D

T_tot = 60.*60.*24.*365. #total mission lifetime, seconds
trials = 999 #number of planets to generate
z = 5 #the critical z score at which a detection is made
t_steps_tot = 100. #total number of times steps that can be spent on a target

top_targets = 20 #how many of the top targets should be printed out

print trials
print scope
print t_steps_tot

file_ = "hygdata_v3.csv" #file contain stellar info

#loads our data and returns the proper p_inside function
#this function gives the probability that a given star has a 
#habitable zone exoplanet inside the field of view for each telescope as a function of x
def load_info():
    x, p_inside_L, p_inside_J, p_inside_W, p_inside_H = pickle.load(open("p_inside_all.p", "rb"))
    probabilities = {"WFIRST":p_inside_W, "JWST":p_inside_J, "LUVOIR":p_inside_L, "HabEx":p_inside_H}
    return x, probabilities[scope]

#finds the fraction of the flux inside the bandwidth
#uses a blackbody approximation
def find_f(T):
    C = 2*math.pi*Math.h*dl*(Math.c**2)/float(Math.steph_boltz*(T**4)*(l**5))
    fac = math.exp(Math.h*Math.c/float(l*T*Math.k)) - 1.
    if fac == 0:
        fac = Math.h*Math.c/float(l*T*Math.k)
    return C/float(fac)


#goes through the stellar file creating a star object for each star 
def get_stars(x, p_inside):
    stars = []
    with open(file_, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        i = 0
        for row in reader:
            if i > 1 and row[16] and row[9] and row[33]:
                star_ = ro.star() #making a new star object
                star_.d = float(row[9])*Math.pc #distance to system
                star_.L = float(row[33])*Math.L_sun #stellar luminosity in W
                root_lum = math.sqrt(float(row[33])) #luminosity in solar luminosities
                star_.x_in = float(star_.d)/float(Math.au*root_lum*1.68) #intrinsice observation parameter
                star_.x = star_.x_in*(l/float(D)) #specific observation parameter
                if star_.d < 30.*Math.pc and star_.x_in < .005e9:
                    ci = float(row[16]) #color index
                    ind = Math.binarySearchClose(star_.x_in/1.0e6, x)
                    star_.p_ = p_inside[ind] #proportion of p(s) that is inside the field of view for this star
                    T_ = 4600.*(1./float(.92*ci + 1.7) + 1./float(.92*ci + .62)) #temperature
                    f_L = find_f(T_) #fraction of luminosity in the bandpass
                    star_.f_L = f_L
                    star_.name = row[4] + " " + row[5] + " " + row[6] + " " + row[15] #some names for the star
                    if star_.p_ > 0: #if the star has >0% chance of having an observable habitable zone exoplanet
                        star_.set_tau(D, l, T, z, C)
                        star_.set_mass()
                        stars.append(star_)
            i += 1
    return stars
    

#simulates trails planets around a given star and integrates until a maximum time has passed or the planet is detected
#returns the 95 percentile integration time and the proportion of planets that were found at all
#both return values here are approximate, with accuracy that is going to increase as root(trials)
#the former variable should be approximately Gaussian distributed, the latter binomial
def find_tau_ex(star_):
    tau, found = [], 0 #a matrix of measured integration times, number of found planets
    for i in range(trials): #planet loop
        star_.add_planet() #create a planet around the star
        tau_, signal, background, z_ = 0, 0, 0, 0 #time integrated, signal photons, background photons, z score
        t_step = star_.planet.T/t_steps_tot #the time steps for the simulation
        while(z_ < z and tau_ < T_tot): #observing loop
            light = star_.get_photons(D, C, l, IWA, OWA, t_step, T) #collect light
            background += light[0]
            signal += light[3] + light[2]
            z_ = (signal - background)/float(math.sqrt(background)) #use cumulative light to calculate z score
            star_.planet.progress(t_step) #progress the planet's orbit
            tau_ += t_step 
        if z_ > z: #check if detected
            tau.append(tau_/T_tot)
            found += 1
        else:
            tau.append(1.0)
    return tau[(len(tau)*19)/20], found/float(trials)

if __name__ == "__main__":
    tau_ex, found = [], [] #estimates of tau exclude for each star (mission lifetimes), estimate of planets that can be found in a mission lifetime for each star
    x, p_inside = load_info()
    stars = get_stars(x, p_inside)
    start = time.time()
    start_ = start
    i = 0
    for star_ in stars: #star loop
        if i%(len(stars)/10) == 0:
            dT = abs(start - time.time())
            print int(((i + 1)/float(len(stars)))*100), " percent complete"
            print abs(start_ - time.time())/(60.*60.), " hours elapsed"
            print dT*abs(trials - i - 1)/(60*60.*(trials/10)), " hours remaining\n"
            start = time.time()
        if star_.p_ > 0:#check to make sure that star is worth looking at
            results = find_tau_ex(star_) 
            star_.tau_exclude = results[0]
            star_.findalbe = results[1]
            tau_ex.append(results[0])
            found.append(results[1])
        i += 1
    stars.sort()
    for i in range(min(top_targets, len(stars))):
        print stars[i].name
    plt.figure(1)
    plt.xlabel("tau to find 95% of findable planets")
    plt.ylabel("counts")
    plt.hist(np.log10(tau_ex), bins = 100)
    plt.show()
    plt.figure(2)
    plt.xlabel("proportion of planets that are findable")
    plt.ylabel("counts")
    plt.hist(found, bins = 100)
    plt.show()
    plt.figure(3)
    plt.xlabel("tau to find 95% of findable planets")
    plt.ylabel("proportion of planets that are findable")
    plt.scatter(np.log10(tau_ex), found)
    plt.show()
    #pickle.dump(stars, open("stars_" + scope + ".p", "wb"))