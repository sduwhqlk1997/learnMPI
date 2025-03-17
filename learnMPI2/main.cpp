#include <iostream>
#include <mpi.h>
using namespace std;

#define MASTER 0
#define MESSAGE 1

void master(int np, int rank)
{
    char message[100];

    sprintf(message,"Broadcast");
    MPI_Bcast(message,100,MPI_CHAR,0,MPI_COMM_WORLD);

    for (int i = 1; i < 10; i++)
    {
        sprintf(message, "Hello %d", i);
        MPI_Send(message, 100, MPI_CHAR, i%(np-1)+1, MESSAGE, MPI_COMM_WORLD);
    } // 发送信息

    for (int i = 1; i < np; i++)
    {
        sprintf(message, "q");
        MPI_Send(message, 100, MPI_CHAR, i, MESSAGE, MPI_COMM_WORLD);
    }
}

void slave(int np, int rank)
{
    char message[100];
    int i;

    MPI_Bcast(message,100,MPI_CHAR,0,MPI_COMM_WORLD);
    printf("Message received bt %d: %s\n",rank,message);

    for (;;)
    {
        MPI_Recv(message, 100, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (message[0] == 'q')
        {
            break;
        }
        printf("Message received bt %d: %s\n", rank, message);
    }
}

int main(int argc, char *argv[])
{
    int np, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // 从MPI_COMM_WORLD中取出np和rank

    if (rank == MASTER)
    {
        master(np, rank);
    }
    else
    {
        slave(np, rank);
    }
    // cout << "testtest" << endl;
    MPI_Finalize();

    return 0;
}