#!/bin/bash

echo "Compiling..."
make
echo "Makefile OK!"

srun --partition=cpar perf stat -e instructions,cycles make runseq CP_CLUSTERS=4
n=4

while [ $n -le 47 ]
do
	echo "running with $n THREADS"
    export OMP_NUM_THREADS=$n
	srun --partition=cpar perf stat -e instructions,cycles make runpar CP_CLUSTERS=4 THREADS=$n

	sleep 5
	n=$((n+4))
done

echo "Script ran flawlessly!"
