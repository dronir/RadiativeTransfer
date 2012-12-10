#!/usr/bin/env python
#coding: utf8

from matplotlib import pyplot as plot
from numpy import loadtxt, degrees, linspace

angles = loadtxt("angles.txt")
intensity = loadtxt("intensity.txt")
plot.hist(angles, 30, color="black")
plot.figure()
plot.bar(linspace(0,180,31), intensity/intensity.sum(), width=180/30)
plot.show()

