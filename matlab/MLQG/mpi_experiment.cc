#include <mpi.h>
#include <iostream>
#include <sstream>
#include <cstdlib>

int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);

    int numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    std::cout << "Running experiment on pid " << pid
              << ", procs = " << numProcs << std::endl;

    std::stringstream command;
    command << "./experiment " << pid << " " << numProcs
            << " > data/experiments/logdir/log" << pid << ".txt";
    std::cout << command.str() << std::endl;
    std::system(command.str().c_str());
    std::cout <<  "./experiment " << pid << " " << numProcs << " finished" << std::endl;
    MPI_Finalize();

    return 0;
}

