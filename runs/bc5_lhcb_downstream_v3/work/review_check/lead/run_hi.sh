#!/bin/bash
cd "$(dirname "$0")"
M=2.2,10.5834,28.527411,62.49
D:/py/python.exe impact.py --variant base --masses $M --nsamp 6000 --ndecay 8 --tag base_hi > hi_base.log 2>&1 &
D:/py/python.exe impact.py --variant rot  --masses $M --nsamp 6000 --ndecay 8 --tag rot_hi  > hi_rot.log 2>&1 &
wait
echo HIDONE
