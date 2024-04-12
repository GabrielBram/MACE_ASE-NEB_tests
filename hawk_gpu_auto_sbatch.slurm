#!/bin/bash --login 
# This submission script was automatically generated using the Calculation Generator at Tue Oct 18 15:28:51 2022. 
# This script is intended for submission to the HAWK_AMD HPC. 
# For any submission script related problems, try fiddling with the settings. 

#SBATCH -n 40
#SBATCH --tasks-per-node=40
#SBATCH -t 0-24:00
#SBATCH --gres=gpu:2
#SBATCH --partition=gpu,gpu_v100
#SBATCH --account=scw1057

# Load the environment
module load python/3.9.2
module load mpi

# module approach to setting up working environment
module use /apps/local/modules/projects
module load scw1057/fhi-aims
module load CUDA/10.2

clush -w $SLURM_NODELIST "sudo /apps/slurm/gpuset_0_shared"

ulimit -s unlimited

# This is introduced as the default PYTHONPATH was changed in March 2022
# Provides backwards compatibility.
# export PYTHONPATH=$PYTHONPATH:/apps/languages/python/3.7.0/el7/AVX512/intel-2018/lib/python3.7/site-packages
. /scratch/c.sacgb4/.python-3.9.2-amd-chemsh/bin/activate
export PYTHONPATH=/home/c.sacgb4/carmm_test/:$PYTHONPATH

export HOME_DIRECTORY=$PWD
#export VERSION=200302
#export VERSION=200511
#export VERSION=210309
#export VERSION=210618
export VERSION=230612
export HOME_DIRECTORY=$PWD
cd $HOME_DIRECTORY


OUTPUT_FILE="test.log"
export output=$OUTPUT_FILE.log
export error=$OUTPUT_FILE.err

time python3 neb_tester.py 1> >(tee $output ) 2> >(tee $error >&2 )
