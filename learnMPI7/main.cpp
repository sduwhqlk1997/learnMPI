#include <iostream>
#include <unistd.h>
#include <mpi.h>
using namespace std;
 int main(int argc, char *argv[]){
    # ifdef DEBUG
    {
        int i=0;
        while (0==i){
            sleep(1);
        }
    }
    #endif
    int numtasks, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (numtasks !=2){
        cout<<"the number of tasks should equal to 2!"<<endl;
        MPI_Finalize();
        exit(0);
    }
    if (rank==0){
        int a =1;
        cout<<"这是主进程"<<endl;
    }
    if(rank==1){
        int a=2;
        cout<<"这是子进程"<<endl;
    }
    MPI_Finalize();

 }
