#!/usr/bin/env python
#coding: utf8

from matplotlib import pyplot as plot
from numpy import loadtxt, degrees, linspace, sin, cos, pi
from sys import argv
from os.path import join

if len(argv) > 1:
    base = argv[1]
else:
    base = "."

angles = loadtxt(join(base, "angles.txt"))
intensity = loadtxt(join(base, "intensity.txt"))
orders = loadtxt(join(base, "orders.txt"))

T = linspace(0,180.0,31)

plot.hist(angles*180/pi, 30, color="black")
plot.figure()
plot.bar(T, intensity, width=180.0/31)
plot.figure()
plot.hist(orders, int(orders.max()))
plot.show()

