#!/usr/bin/env python

import sys
import random

m = int(sys.argv[1])
n = int(sys.argv[2])

i = 0

f = open("matrix.txt","w")
v = open("vector.txt","w")

f.write(str(n) + "\n")
v.write(str(n) + "\n")

while (i < n):
    j = 0;
    v.write(str(random.randint(0,1000)*random.randint(0,m)) + "\t")
    while (j < n) :
        if (i-m/2 <= j <= i+m/2) :
            f.write(str(random.randint(0,1000)*random.randint(0,m)) + "\t")
        else:
            f.write("0 \t")
        j += 1
    f.write("\n")
    i += 1

