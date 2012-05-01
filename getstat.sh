#!/bin/bash

rm ./output.txt
rm ./toutput.txt
i=0

# $1 - number of test
# $2 - width
# $3 - dimension
# $4 - type of matrix
# $5 - pivtol

while [ $i -lt $1 ]
do
  ./generator $2 $3 $4
  ./sparce3 $5
  i=$[$i+1]
  echo $i
done

python ./stat.py $3
