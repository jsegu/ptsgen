#!/bin/bash

for func in ramp cos epica grip odp1012
do
	python2 ptsgen.py $func -120000 0 0 -10 -o $func.nc
done

./plot.py
