#!/bin/bash

file="RiemannExactResults/Test3/Test3_100cells.dat"

for i in {2..5}
do
	echo "plot '$file' using 1:$i with linespoints" | gnuplot -persist
done	
