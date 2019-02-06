#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=mcsim_1          ## sets the job name
#SBATCH --time=4-00:00:00           ## sets to 4 day the job's runtime wall-clock limit.
#SBATCH --nodes=8		            ## assigns 4 cores for execution.
#SBATCH --ntasks-per-node=1         ## Request 1 core per node.
#SBATCH --mem=64M	                ## sets to 32MB as the per process enforceable memory limit.
#SBATCH --output=apap58.1107.%j	        ## directs the job's standard output (%J is jobid)
#SBATCH --mail-type=END 
#SBATCH --mail-user=nhsieh@cvm.tamu.edu    #Send all emails to email_address
#SBATCH --account=122803364714     ## Change account

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/nhsieh/lib
export LD_LIBRARY_PATH
source ~/.bash_profile

## For sensitive parameters
time ./mcsim.apap.pbpk_v2 apap58.c01.in &
time ./mcsim.apap.pbpk_v2 apap58.c02.in & 
time ./mcsim.apap.pbpk_v2 apap58.c03.in & 
time ./mcsim.apap.pbpk_v2 apap58.c04.in &
time ./mcsim.apap.pbpk_v2 apap58.c05.in &
time ./mcsim.apap.pbpk_v2 apap58.c06.in & 
time ./mcsim.apap.pbpk_v2 apap58.c07.in & 
time ./mcsim.apap.pbpk_v2 apap58.c08.in
wait