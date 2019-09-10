#module load slurm anaconda3
source activate rnaseq

# Display all variables set by slurm
env | grep "^SLURM" | sort

unset DISPLAY
