#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

f = open("output.txt","r")

a = []
b = [0,0,0,0,0]
c = ["MAX_MULT","MAX_SUM","NORM_EXACT","NORM_MULT","NORM_SUM"]
q = 0;
for line in f.xreadlines():
    a = line.split(" ");
    for i in range(0,5):
        b[i] += float((int(a[i+1])-int(a[0])))/int(a[0])
    q +=1;

for i in range(0,5):
    b[i] = b[i]*100/q;
    print c[i] , "\t %0.1f" % b[i], "%"