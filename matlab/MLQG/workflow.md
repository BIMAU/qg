# Workflow: how to get a parallel experiment running on Cartesius

1. Compile `mpi_experiment.cc` (probably just once):
   `mpicxx -g -Wall mpi_experiment.cc -o mpi_experiment -lmpi`

2. Edit and compile `experiment.m`:
   `./compile_mfile.sh experiment.m`
   
3. Edit and submit `submit_mpi_experiment.sh`:
   `sbatch submit_mpi_experiment.sh`
