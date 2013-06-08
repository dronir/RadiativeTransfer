#!/usr/bin/env python
#coding: utf8

from pylab import *
from sys import argv

data = loadtxt(argv[1], delimiter=",")

alpha = data[:,0]
first = data[:,1]
total = data[:,2]
analytic = data[:,3]

plot(alpha, first, color="red", linewidth=2, label="first-order")
plot(alpha, total, color="blue", linewidth=2, label="total")
plot(alpha, analytic, "--", color="black", linewidth=2, label="analytical")
ylim(0,total.max())
xlim(0,pi)
legend()
show()
