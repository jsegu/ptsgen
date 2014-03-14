#!/bin/bash

for func in ramp cos lapaz21p odp1012 odp1020 domefuji epica vostok grip
do
	python2 ptsgen.py $func -120000 0 0 -10 -o $func.nc
done

./plot.py
