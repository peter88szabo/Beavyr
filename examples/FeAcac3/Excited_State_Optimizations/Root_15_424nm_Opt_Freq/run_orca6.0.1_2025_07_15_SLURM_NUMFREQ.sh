#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=72:00:00
#SBATCH -A lp_harvey_first
#SBATCH -p batch_sapphirerapids
#SBATCH -J R15NFreq


FILE="Exc_root_OnlyFreq_from_optimized"
GEOM="Exc_root_opt_continue.xyz"
WAVEFUNC="Exc_root_opt_continue.gbw"
ORCA_LOCAL_PLACE="/data/leuven/350/vsc35002/Programs/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2"
WORKDIR=$VSC_SCRATCH_NODE
HOME=`pwd`

echo "Home: " $HOME
echo "Workdir: " $WORKDIR

#Either use the OpenMPI from the module library:
module load OpenMPI/4.1.5-GCC-12.3.0
#or this locally installed mpi_conda_env
#cp -r $VSC_DATA/mpi_conda_env $WORKDIR/
#export PATH=$PATH:$WORKDIR/mpi_conda_env/bin
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORKDIR/mpi_conda_env/lib

#optinally copy the whole orca dir to the scratch disk together witht the mpi_conda_env
#cp -rf ${ORCA_LOCAL_PLACE} ${WORKDIR}
#orca_run=orca_6_0_1_linux_x86-64_shared_openmpi416_avx2/orca
orca_run=${ORCA_LOCAL_PLACE}/orca

#-------------------------------------------------
backup_files() {
    while true; do
        sleep 300  # Wait 5 minutes
        cp -u *.out ${HOME}
    done
}
#-------------------------------------------------

cp ${FILE}.inp ${WAVEFUNC} ${GEOM} $WORKDIR

cd $WORKDIR

# Start backup process in the background and store its ID
backup_files &
BACKUP_PID=$!

$orca_run ${FILE}.inp > ${FILE}.out

rsync -av . $HOME --exclude=*.tmp --exclude=*.inp --exclude='mpi_conda_env/' 

kill $BACKUP_PID
