#!/bin/bash
#$ -q x19.q
#$ -pe x40 40
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -N FNAME.8c
#

module load intel/2020.2.254
module load intelmpi/2020.2.254
module load intelMKL/2022.0.2
source /home/opt/settings/2017.4/intel-mpi.sh
module load qe/7.2

I_MPI_PIN=1
OMP_NUM_THREADS=1
PW_COMMAND=/home/krojas/share/software/qe/qe-7.2/bin
MPI_COMMAND="mpirun -n $NSLOTS"

###ln -fs ~/QE/qe-6.3/bin/vdW_kernel_table

prefix1=8c
$MPI_COMMAND ${PW_COMMAND}/pw.x --npool $NSLOTS < $prefix1.pwi > $prefix1.pwo
