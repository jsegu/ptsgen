#!/bin/bash
# Copyright (c) 2014--2018, Julien Seguinot <seguinot@vaw.baug.ethz.ch>
# GNU General Public License v3.0+ (https://www.gnu.org/licenses/gpl-3.0.txt)

for func in ramp cos lapaz21p odp1012 odp1020 domefuji epica vostok ngrip grip guliya md012444
do
	python ptsgen.py $func -120000 0 0 -10 -o $func.nc
done

for expt in {0..3}
do
	python ptsgen.py grip -120000 0 0 -10 -s $((10**$expt)) -o grip-1e$expt.nc
done

./plot.py
