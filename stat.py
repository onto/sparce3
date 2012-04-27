#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

f = open("output.txt","r")
tf = open("toutput.txt","r")

a = []
b = [0,0,0,0,0,0,0,0]
c = ["MAX_EXACT","MAX_SUM","MAX_MULT","MAX_MULT2","NORM_EXACT","NORM_SUM","NORM_MULT","NORM_MULT2"]
#c = ["MAX_SUM"]
q = [0,0,0,0,0,0,0,0]
z = [0,0,0,0,0,0,0,0]
t = [0,0,0,0,0,0,0,0]
tc = [0,0,0,0,0,0,0,0]
nz = 0;
n = int(sys.argv[1])
j = 0;


for line in f.xreadlines():
    a = line.split(" ");
    nz += float(a[0])/(n*n)
    for i in range(0,8):
        if (int(a[i+1]) != 0):
            b[i] += (float(a[i+1])-float(a[0]))/float(a[0])
            q[i] += 1;
            z[i] += float(a[i+1])/(n*n)
    j += 1

for line in tf.xreadlines():
    a = line.split(" ");
    for i in range(0,8):
        if (float(a[i]) != -1):
            t[i] += float(a[i])
            tc[i] += 1

nz = nz*100/j
print "NZ: %0.2f" % nz , "%"

for i in range(0,8):
    if (q[i] != 0) :
        b[i] = b[i]*100/q[i];
        z[i] = z[i]*100/q[i];
        if (tc[i] != 0):
            t[i] = t[i]/tc[i];
        print c[i] , "\t %0.1f" % b[i], "%", "\t %0.1f" % z[i], "%", "\t %0.3f" % t[i]
