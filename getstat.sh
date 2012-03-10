#!/bin/bash

rm ./output.txt
i=0
while [ $i -lt $1 ]
do
  python ./generator.py $2 $3
  ./sparce3
  i=$[$i+1]
done

python ./stat.py
