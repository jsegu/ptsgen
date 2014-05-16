#!/bin/bash

for func in ramp cos lapaz21p odp1012 odp1020 domefuji epica vostok ngrip grip guliya
do
	python2 ptsgen.py $func -120000 0 0 -10 -o $func.nc
done

./plot.py
