# -*- coding: utf-8 -*-
"""
Created on Mon Mar 06 12:22:22 2017

@author: aeber
"""
#contains the orbit class which generates random orbital parameters
import random
import Math
import math
import numpy as np

#randomly picks and returns an eccentricity and semi major axis
#given minimum semi major axis (a_min), maximum semi major axis (a_max) and maximum eccentricity (e_max)
def get_ea(a_min, a_max, e_max):
    looking = True
    while(looking):
        e_ = random.uniform(0, e_max)
        a_ = random.uniform(a_min, a_max)
        if a_*(1. - e_) > a_min and a_*(1. + e_) < a_max:
            return e_, a_
            looking = False

#returns the 3-space seperation pdf values given a semi major axis (a), eccentricity (e), and 3-space seperation (r)
def p_r(a, e, r):
    C = 1./(a**(2./3.)*math.pi)
    num = math.sqrt((r*e)**2 + 2.*r*a*(1. - e**2) - r**2)
    denom = math.sqrt(((r*e)**2 -(a*(1. - e**2) - r)**2)*(2./r + 1./a))
    return float(C*num)/denom

#given original coordinates (x,y,z), rotation angle (th), axis = [x, y, z]
#returns new coordinate values (x_, y_, z_)
def rotate(x,y,z, th, axis = 2):
    cos, sin = math.cos, math.sin
    x_, y_, z_ = 0, 0, 0
    if axis == 1:
        x_ = cos(th)*x + sin(th)*z
        y_ = y
        z_ = -1*sin(th)*x + cos(th)*z
    elif axis == 2:
        x_ = cos(th)*x - sin(th)*y
        y_ = sin(th)*x + cos(th)*y
        z_ = z
    else:
        x_ = x
        y_ = -1*sin(th)*z + cos(th)*y
        z_ = cos(th)*z + sin(th)*y   
    return x_, y_, z_ 
    

#randomly picks and returns a 3-space seperation given semi major axis (a) and eccentricitry (e)
def get_r(a, e):
    looking = True
    while(looking):
        max_ = 1./(a*2.*math.pi*e)
        y = random.uniform(0, 100*max_)
        x = random.uniform(a*(1. - e),a*(1. + e))
        if y < p_r(a, e, x):
            looking = False
            return x

#returns the distance given x, y
def distance(x, y):
    return math.sqrt(x**2 + y**2)

#converts a 3-space sepration value to a true anomaly value
def r2nu(r, a, e):
    arg = a*(1. - e**2)
    arg /= float(r*e)
    arg -= 1./e
    return math.acos(arg)

#converts a true anomaly value to a 3-space seperation value  
def nu2r(nu, a, e):
    return (a*(1. - e**2))/float(1. + e*math.cos(nu))

    
#the orbit class
class orbit(object):
    
    #pass luminosity in W, mass in kg
    def __init__(self, lum, M, rand = True):
        self.lum_star = lum #luminosity of host star in Watts
        self.M = M #mass of host star in kg
        self.a = 0 #semi major axis, m
        self.e = 0 #eccentiricty
        self.T = 0 #period, s
        self.thetas = [] #randomly generated rotation angles
        self.axes = [] #randomly generated rotation axes
        self.nu = 0 #true anomaly
        self.r = 0 #three space seperation
        self.R = Math.r_e #radius, m
        #coordinates in plane of orbit
        self.x = 0 
        self.y = 0
        #coordinates in plane of observation
        self.x_ = 0
        self.y_ = 0
        self.z_ = 0
        self.a_max = 1.68*math.sqrt(lum/float(Math.L_sun))*Math.au #maximum  possible semi major axis
        self.a_min = math.sqrt(lum/float(Math.L_sun))*Math.au #minimum possible semi major axis
        self.e_max = (self.a_max - self.a_min)/float(self.a_max + self.a_min) #maximum possible eccentricity
        if rand: #randomly orient self
            self.randomize()
    
    #gives the orbit random parameters
    def randomize(self):
        self.e, self.a = get_ea(self.a_min, self.a_max, self.e_max)
        self.T = 2.*math.pi*math.sqrt((self.a**3)/float(Math.G*self.M))
        for i in range(10):
            theta_ = random.uniform(0, 2.*math.pi)
            axis_ = random.randint(0,2)
            self.thetas.append(theta_)
            self.axes.append(axis_)
        self.r = get_r(self.a, self.e)
        self.nu = r2nu(self.r, self.a, self.e)
        self.set_xy()
   
    #evolves the orbit, (i.e. this is progress the verb not the noun) 
    def progress(self, dt, dth = math.pi/100.):
        dA_tot = (dt/float(self.T))*math.pi*self.a*self.a*(1. - self.e**2)
        dA = 0
        while(dA < dA_tot):
            dA += (.5*(self.r)**2)*dth
            self.nu += dth
            self.r = nu2r(self.nu, self.a, self.e)
        self.set_xy()
    
        
    #sets the coordinate variables ot the correct values
    def set_xy(self):
        self.x = self.r*math.cos(self.nu)
        self.y = self.r*math.sin(self.nu)
        X, Y, Z = self.r*math.cos(self.nu), self.r*math.sin(self.nu), 0
        for i in range(len(self.thetas)):
            X, Y, Z = rotate(X, Y, Z, self.thetas[i], self.axes[i])
        self.x_, self.y_, self.z_ = X, Y, Z


#object containing information on star            
class star(object):
    
    def __init__(self, L = Math.L_sun, d = Math.pc, x_in = 1., x = 1., p_ = 0):
        self.name = "" #name of the star
        self.L = L #luminsoty in watts
        self.d = d #distance to system in m
        self.x_in = x_in #intrinsic observation parameter
        self.x = x #specific observation parameter
        self.planet = None #planet orbit object around this star
        self.tau = 0 #expected necessary observation time (this is just a proxy), s
        self.tau_exclude = 0 #measured time necessary to exclude 95% HZExs
        self.findalbe = 0 #measured proportion of exoplanets that are findable
        self.M = 0 #mass in kg
        self.p_ = p_ #probability of being inside field of view
        self.f_L = 0 #fraction of luminosity inside bandwidth
        self.t_exclude = 0 #time required to exclude 95% of HZ orbits, if 0 then it is impossible in a mission lifetime
    
    #compares two stars, precedence is taken by star with shorter tau_exclude
    def __lt__(self, other):
        use_findable = max(1e-5, self.findalbe)
        other_find = max(1e-5, other.findalbe)
        if self.tau_exclude/use_findable < other.tau_exclude/other_find:
            return True
        elif self.t_exclude/use_findable > other.tau_exclude/other_find:
            return False
        if self.tau < other.tau:
            return True
        return False
    
    #find the expectation value for the observation time
    #given aperture diameter, observation wavelength, transmission, z score, contrast
    def set_tau(self, D, l, T, z, C):
        A_phi = .15 #average phase multiplied by average albedo
        A = (Math.au/float(Math.r_e))**4
        A *= (z*self.d/float(D*self.p_*Math.L_sun*A_phi))**2
        A *= (C*Math.h*Math.c*self.L)/float(T*l*self.f_L)
        self.tau = A
    
    #sets the mass of the star using an approximate mass-luminosity relation
    def set_mass(self):
        l_ = self.L/float(Math.L_sun)
        if l_ < .033:
            self.M = Math.M_sun*((l_/.23)**(1./2.3))
        elif l_ < 16.:
            self.M = Math.M_sun*((l_)**(1./4.))
        elif l_ < 53700:
            self.M = Math.M_sun*((l_/1.5)**(1./3.5))
        else:
            self.M = Math.M_sun*l_/3200.
            
    #generates an oribt around this star
    def add_planet(self):
        self.planet = orbit(self.L, self.M)
    
    #returns background, signal (expectation values), and background, signal randomly chosen
    #given aperture diameter, contrast, wavelength, inner working angle in l/D
    #Outer working angle l/D, integration time, transmission
    def get_photons(self, D, C, l, IWA, OWA, dt, T):
        L_ = (self.f_L*self.L*dt)/float(Math.h*Math.c/float(l)) #photons from star
        back = T*C*L_*(D**2)/float(4*self.d**2)
        s_ = math.sqrt(self.planet.x_**2 + self.planet.y_**2)/Math.au
        field = (s_/float(self.d/Math.au) > IWA*l/float(D)) and (s_/float(self.d/Math.au) < OWA*l/float(D))
        try:
            phase = (math.acos(s_/float(self.planet.r/Math.au)) + math.pi/2.)/float(math.pi)
        except ValueError:
            phase = .5
        if self.planet.z_ > 0:
            try:
                phase = (math.pi/2. - math.acos(s_/Math.au/float(self.planet.r/Math.au)))/float(math.pi)
            except ValueError:
                phase = .5
        signal = (field)*phase*T*L_*((self.planet.R*D)**2)/(float(16*self.planet.r*self.d)**2)
        signal_ = signal
        back_ = back
        if signal_ < 1e4 and signal_ > 0:
            signal_ = np.random.poisson(signal_)
        if back_ < 1e4:
            back_ = np.random.poisson(back_)
        return back, signal, back_, signal_
            
        