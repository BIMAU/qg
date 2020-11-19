# Workflow: how to get a parallel experiment running on Cartesius

1. Compile `mpi_experiment.cc` (probably just once):
   `mpicxx -g -Wall run_trans_mpi.cc -o run_trans_mpi -lmpi`

2. Edit and compile `experiment.m`:
   `./compile_mfile.sh experiment.m`
   
3. Edit and submit `submit_mpi_experiment.sh`:
   `sbatch submit_mpi_experiment.sh`
