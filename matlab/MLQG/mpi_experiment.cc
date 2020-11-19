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

    std::cout << "Running transient on pid " << pid
              << ", numProcs = " << numProcs << std::endl;

    std::stringstream create_dir;
    create_dir << "mkdir -p data/experiments/logdir";
    std::system(create_dir.str().c_str());
    
    std::stringstream command;
    command << "./experiment " << pid << " " << numProcs
            << " > data/experiments/logdir/log" << pid << ".txt";
    std::cout << command.str() << std::endl;
    std::system(command.str().c_str());

    MPI_Finalize();

    return 0;
}

