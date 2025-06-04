#!/bin/bash
#$ -q x20.q
#$ -pe x52 52
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -N z6-2-1
#$ -t 1-30:1
#
source /home/opt/settings/2017.4/intel-compiler.sh
source /home/opt/settings/2017.4/intel-mpi.sh

export PYTHONPATH=/home/bctan/GPAW/src/gofeeNEW:$PYTHONPATH
export PATH=$HOME/GPAW/src/gofeeNEW:$PATH
export PATH=$HOME/GPAW/src/gofeeNEW/statistics_tools/survival_statistics:$PATH

I_MPI_PIN=1
OMP_NUM_THREADS=1
I_MPI_ADJUST_ALLGATHERV=2

id=${SGE_TASK_ID}
id="$(($id-1))"

#EDIT
mkdir runs0
dir="runs0/run$id"
pyr=run_search.py

mkdir -p $dir
cp $pyr $dir
cd $dir

mpirun -np $NSLOTS python $pyr > search.log
