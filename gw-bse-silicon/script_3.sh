#!/bin/bash -l
#SBATCH -t 00:30:00
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -c 1

# Edit these lines as needed
BGW_BIN="$HOME/bgw/BerkeleyGW-4.0/bin"
EPSILON="$BGW_BIN/epsilon.real.x"
SIGMA="$BGW_BIN/sigma.real.x"
KERNEL="$BGW_BIN/kernel.real.x"
ABSORPTION="$BGW_BIN/absorption.real.x"
export OMP_NUM_THREADS=1
# put argument for number of procs here too if needed, e.g. -n 8
# MPIRUN="srun -n 32 -c 1"
MPIRUN="mpirun -np 24"

#
cd ./07-epsilon
$MPIRUN $EPSILON &> ./epsilon.out
echo "07-epsilon DONE"
cd ..
#
cd ./08-sigma
$MPIRUN $SIGMA &> ./sigma.out
echo "08-sigma DONE"
cd ..
#
cd ./09-kernel
$MPIRUN $KERNEL &> ./kernel.out
echo "09-kernel DONE"
cd ..
#
cd ./10-absorption
ln -sf ../08-sigma/eqp1.dat eqp_co.dat
$MPIRUN $ABSORPTION &> ./absorption.out
echo "10-absorption DONE"
cd ..
