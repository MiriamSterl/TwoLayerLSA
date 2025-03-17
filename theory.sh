#!/bin/bash -l
#
#SBATCH -J theory
#SBATCH -p normal
#SBATCH -t 120:00:00
#SBATCH -o theory.%j.o
#SBATCH -e theory.%j.e

cd ...

conda activate base

python3 theory_nonzonal.py "topo_nonzonal_U_45deg_nofric_fine" "0" "45" "1"
python3 theory_nonzonal.py "topo_nonzonal_U_300deg_nofric_fine" "0" "300" "1"

python3 theory_nonzonal.py "topo_nonzonal_U_45deg_mu_weak_fine" "100" "45" "1"
python3 theory_nonzonal.py "topo_nonzonal_U_300deg_mu_weak_fine" "100" "300" "1"

python3 theory_nonzonal.py "topo_nonzonal_U_45deg_mu_intm_fine" "10" "45" "1"
python3 theory_nonzonal.py "topo_nonzonal_U_300deg_mu_intm_fine" "10" "300" "1"