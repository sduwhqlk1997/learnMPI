#include <iostream>
#include <mpi.h>
#define NPROCS 8

int main(int argc, char *argv[])
{
    int rank, new_rank, sendbuf, recvbuf, numtasks,
        ranks1[4] = {0, 2, 4, 6}, ranks2[4] = {1, 3, 5, 7};
    MPI_Group orig_group, new_group; // required variables
    MPI_Comm new_comm;               // required variable

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    if (numtasks != NPROCS)
    {
        printf("Must specify MP_PROCS= %d. Terminating.\n", NPROCS);
        MPI_Finalize();
        exit(0);
    }

    sendbuf = rank;

    // extract the original group handle
    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

    // divide tasks into two distinct groups based upon rank
    if (rank % 2 == 0)
    {
        MPI_Group_incl(orig_group, NPROCS / 2, ranks1, &new_group);
    }
    else
    {
        MPI_Group_incl(orig_group, NPROCS / 2, ranks2, &new_group);
    }

    // create new new communicator and the perform collective communications
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    if (new_comm != MPI_COMM_NULL)
    {
        MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, new_comm);

        // get rank in new group
        MPI_Group_rank(new_group, &new_rank);
        printf("rank= %d newrank= %d recvbuf= %d\n", rank, new_rank, recvbuf);
    }

    MPI_Finalize();
}