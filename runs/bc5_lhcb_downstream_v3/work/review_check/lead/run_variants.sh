#!/bin/bash
cd "$(dirname "$0")"
M=0.5,1.0,1.5,2.0,2.2,3.0,4.34457,5.082255,7.5,10.5834,15.234347,28.527411,45.66582,62.49
for v in rot libmass kLdead multi4pi gamma bxi bincl bsext; do
  D:/py/python.exe impact.py --variant $v --masses $M --nsamp 1200 --ndecay 8 > var_$v.log 2>&1 &
done
wait
echo ALLDONE
