#!/bin/bash

rm ./output.txt
rm ./toutput.txt
i=0
while [ $i -lt $1 ]
do
  python ./generator.py $2 $3
  ./sparce3 $4
  i=$[$i+1]
  echo $i
done

python ./stat.py $3
