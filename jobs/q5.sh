#!/usr/local_rwth/bin/zsh
 
#SBATCH -c 24

#SBATCH --time=1:00:00

# name the job
#SBATCH --job-name=HPMC-Exam
 
# declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.txt
##SBATCH --account=lect0070
#SBATCH --gres=gpu:volta:1
 
### beginning of executable commands
file="../data/q5.csv"
# ../main.out --parsable -r 5 10000 >> ${file}


