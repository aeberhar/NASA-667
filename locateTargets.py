# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 09:01:45 2016

@author: aeber
"""

#converts a data table on exoplanets into a pickle file with ToI objects

direct = "C:/Users/aeber/Desktop/Python/NASA"
import sys #module used to add a path to our script so we can add the following modules
sys.path.insert(0, direct)
import csv #a reader for comma seperated values
import candidate1_1 as candidate #the candidate class
import pickle #pickle module, lets me save and output data
import Math

#give index of each parameter, next two indices are +/- error respectively
mIndex = 21 #mass of planet
TIndex = 53 #index for effective temperature
rIndex = 61 #index of steller radius
OIndex = 46 #optical magnitude index
DIndex = 42 #index of distance to star
LIndex = 202 #steller luminosity index
aIndex = 9 #index of semimajor axis
indices_ = [mIndex, TIndex, rIndex, OIndex, DIndex, LIndex, aIndex]

#factors to convert each value given in data to mks units
factors_ = [Math.m_jupiter, 1.0, Math.r_solar, 0, Math.pc, Math.L_solar]


#returns argument x as a float, if typeError then returns 0
def float_(x):
    try:
        return float(x)
    except TypeError:
        return 0.


#outputs the data to a pickle file with given data
def output(data, name):
    ofile = open(name, 'wb')
    pickle.dump
    ofile.close()


#adds a candidate object using the data in row to given planet matrix (planets),
#uses logical cut functions cutFunc to discriminate against objects
#given list of unit conversion factors and attribute indices in the given row
def addCandidate(row, planets, cutFunc, factors, indices):
    o = candidate.candidate(float_(row[0]))
    addFirstAttr(row, o, factors, indices)
    addSecondAttr(row, o, factors, indices)
    o.remainder()
    for i in range(len(cutFunc)):
        o.p.append(cutFunc[i](o))
    o_ = candidate.ToI()
    o_.copy(o)
    planets.append(o_)


#adds mass, effective steller tempetature, steller radius, and steller magnitude attributes
def addFirstAttr(row, o, f, i):
    addAttr(o.m, f[0], row, i[0])
    addAttr(o.T_eff, f[1], row, i[1])
    addAttr(o.r_star, f[2], row, i[2])
    o.mag[0] = float_(row[i[3]])
    o.mag[1] = float_(row[i[3]+1]) 
    

#adds distance, luminosity, and orbital radius attributes
def addSecondAttr(row, o, f, i):
    addAttr(o.D, f[4], row, i[4])
    try:
        o.L[0] = f[5]*(10**float(i[5]))
    except ValueError:
        pass
    try:
        o.L[1] = f[5]*(10**float(i[5]+1))
    except ValueError:
        pass
    o.getL()
    o.getD()
    o.getL()
    addAttr(o.a, f[6], row, i[6])
    o.geta()
    

#looks to see if there is sufficient data to add a given arbitrary attribude to
#our planet object, converts data point row[i] to mks units using given factor
def addAttr(attr, factor, row, i):
    attr[0] = factor*float_(row[i])
    attr[1] = factor*(abs(float_(row[i + 2])) + float_(row[i + 1]))/2.0


#the main file for this script, takes a data fileName (infileName);
#output file name (ofileName), row of infile on which the data starts (startRow)
#cutFunc is a list of functions that take candidate objects and return p values
#factors is an array of factors that convert the infile params to mks units
#indices has the indices of each parameter in each data row
def main(infileName, ofileName, startRow = 0, cutFunc = [], factors = factors_, indices = indices_):
    i = 0
    planets = []
    with open(infileName, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if i >= startRow:
                addCandidate(row, planets, cutFunc, factors, indices)
            i += 1
    output(planets, ofileName)