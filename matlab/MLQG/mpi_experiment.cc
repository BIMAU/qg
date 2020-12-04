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

    std::cout << "Running main on pid " << pid
              << ", procs = " << numProcs << std::endl;

    std::string threads = "1";
    if (argc == 2)
        threads = argv[1];

    std::stringstream command;
    command << "./main " <<  pid << " " << numProcs << " " << threads;
    // << " > data/experiments/logdir/log" << pid << ".txt";s
    std::cout << command.str() << std::endl;
    std::system(command.str().c_str());
    std::cout <<  "./main " << pid << " " << numProcs << " finished" << std::endl;
    MPI_Finalize();

    return 0;
}
