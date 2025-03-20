#include <iostream>
#include <mpi.h>
using namespace std;
#define NRA 100
#define NCA 50
#define NRB 50
#define NCB 30
#define FROM_MASTER 1
#define FROM_WORKER 2
#define MASTER 0
int main(int argc, char *argv[]){
    double A[NRA][NCA];
    double B[NRB][NCB];
    double C[NRA][NCB];

    int numtasks, taskid, offset, tag,rows;
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

    /* Master task */
    if (taskid==MASTER){
        int SendNR = NRA/numtasks; // 传递的行数
        int remainderNR = NRA % numtasks; // 行的余数
        offset = 0;
        tag=FROM_MASTER;
        for(int dest=1;dest<numtasks;dest++){ // 向各子进程传递信息
            rows = dest < remainderNR ? SendNR+1 : SendNR;
            MPI_Send(&A[offset][0],rows*NCA,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
            MPI_Send(&offset,1,MPI_INT,dest,tag,MPI_COMM_WORLD);
            MPI_Send(&rows,1,MPI_INT,dest,tag,MPI_COMM_WORLD);
            offset+=rows;
        }

        tag = FROM_WORKER;
        for(int dest=1;dest<numtasks;dest++){ // 接收子进程的计算结果
            MPI_Recv(&offset, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&C[offset][0], rows*NCB, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &status);
        }
        /* 展示计算结果 */
        for(int i=0;i<NRA;i++){
            for(int j=0;j<NCB;j++){
                cout<<C[i][j]<<"\t";
            }
            cout<<endl;
        }
    }
    
    /* 子进程 */
    if(taskid > MASTER){
        /* 接收主进程的信息 */
        tag=FROM_MASTER;
        
    }
}