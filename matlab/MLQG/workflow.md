# Workflow: how to get a parallel experiment running on Cartesius

1. Edit `experiment.m` and other relevant `.m`-files:
   
2. `./compile_mfile.sh main.m`

3. Compile `mpi_experiment.cc` (probably just once):
   `mpicxx -g -Wall mpi_experiment.cc -o mpi_experiment -lmpi`
   
4. Edit `submit_mpi_experiment.sh`:

5. `sbatch submit_mpi_experiment.sh`
