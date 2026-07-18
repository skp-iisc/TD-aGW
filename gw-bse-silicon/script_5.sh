#!/bin/bash -l
#SBATCH -t 00:05:00
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1

# Note: you can actually run this script in serial if you compiled BerkeleyGW in serial.

# Edit these lines as needed
# module load python
BGW_BIN="$HOME/bgw/BerkeleyGW-4.0/bin"
INTEQP="$BGW_BIN/inteqp.real.x"
export OMP_NUM_THREADS=1
# MPIRUN="srun -n 1 -c 1"
MPIRUN="mpirun -np 20"

# Plot interpolated bandstruture using the inteqp code.
# This step works with both ESPRESSO and PARATEC as starting mean-field codes.
cd 12-inteqp
ln -sf ../08-sigma/eqp1.dat eqp_co.dat
$MPIRUN $INTEQP &> inteqp.out
echo "12-inteqp DONE"
# ./plot_inteqp.py
python plot_inteqp.py
cd ../

cd 13-plots

python inteqp_bands.py
echo "band structure plotted."

python absorption.py
echo "absorption (eh/no eh) spectrum plotted."

cd ..
