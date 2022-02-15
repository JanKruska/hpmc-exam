#!/usr/local_rwth/bin/zsh
 
#SBATCH -c 48

#SBATCH --time=1:00:00

# name the job
#SBATCH --job-name=HPMC-Exam
 
# declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.txt
#SBATCH --account=lect0070
 
### beginning of executable commands
## Set env variables
export OMP_NUM_THREADS=1
file="../data/q3.csv"
../main.out -r 5 --parsable 1000 10000 500 > ${file}
