#!/usr/local_rwth/bin/zsh
 
#SBATCH -c 48

#SBATCH --time=1:00:00

# name the job
#SBATCH --job-name=HPMC-Exam
 
# declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.txt
#SBATCH --account=lect0070
 
### beginning of executable commands
echo "n,m,time,GFLOPS,GFLOPS/core" > q4.csv
for n in 1 2 4 8 16 24 32 48
do
  export OMP_NUM_THREADS=${n}
  echo -n "${n}," >> q4.csv
  ./main.out --parsable --no-header -r 5 10000 >> q4.csv
done

