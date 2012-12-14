#!/usr/bin/env python
#coding: utf8

from matplotlib import pyplot as plot
from numpy import loadtxt, degrees, linspace, sin, cos, pi

angles = loadtxt("angles.txt")
intensity = loadtxt("intensity.txt")
orders = loadtxt("orders.txt")

T = linspace(0,180.0,31)

plot.hist(angles*180/pi, 30, color="black")
plot.figure()
plot.bar(T, intensity, width=180.0/31)
plot.figure()
plot.hist(orders, int(orders.max()))
plot.show()

