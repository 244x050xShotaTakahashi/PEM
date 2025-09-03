#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr20001a
#SBATCH -t 72:00:00
#SBATCH --rsc p=2:t=2:c=4:m=10000M
#SBATCH -o %x.%j.out

#============ Shell Script ============
date
echo ...starting job...
/opt/system/app/intelpython/2024.2.0/bin/python animate_pem.py

date