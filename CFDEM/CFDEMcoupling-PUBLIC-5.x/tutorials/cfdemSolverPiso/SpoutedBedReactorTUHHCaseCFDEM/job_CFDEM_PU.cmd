#!/bin/bash
# parallel job using 16 cores. and runs for 4 hours (max)
#SBATCH -N 1                 # Request that a minimum of nodes be allocated to this job
#SBATCH --ntasks-per-node=16 # Request the maximum ntasks be invoked on each node
#SBATCH -t 360:00:00           # Set a limit on the total run time of the job allocation in "hours:minutes:seconds"
#SBATCH --mem=64000            # Specify the real memory required per node in MegaBytes
# sends mail when process begins, and
# when it ends. Make sure you define your email
# address.
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=yourNetID@princeton.edu
#SBATCH -J TUHHBox

#Run case
#source /home/hema/mostafas/OpenFOAM/OpenFOAM-2.2.x/etc/bashrc
./parCFDDEMrun.sh

