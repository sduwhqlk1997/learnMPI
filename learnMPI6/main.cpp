#include <iostream>
#include <mpi.h>
using namespace std;
#define NRA 20
#define NCA 10
#define NRB 10
#define NCB 10
#define FROM_MASTER 1
#define FROM_WORKER 2
#define MASTER 0
int main(int argc, char *argv[])
{
    double A[NRA][NCA];
    double B[NRB][NCB];
    double C[NRA][NCB];

    int numtasks, taskid, offset, tag, rows;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    /* Master task */
    if (taskid == MASTER)
    {
        cout << "MPI_MM has started with " << numtasks << " tasks" << endl;
        numtasks -= 1;
        cout << "Initializing arrays..." << endl;
        for (int i = 0; i < NRA; i++)
        {
            for (int j = 0; j < NCA; j++)
            {
                A[i][j] = i + j;
            }
        } // 初始化A
        for (int i = 0; i < NRB; i++)
        {
            for (int j = 0; j < NCB; j++)
            {
                B[i][j] = i * j;
            }
        } // 初始化B

        int SendNR = NRA / numtasks;      // 传递的行数
        int remainderNR = NRA % numtasks; // 行的余数
        offset = 0;
        tag = FROM_MASTER;
        for (int dest = 1; dest <= numtasks; dest++)
        { // 向各子进程传递信息
            rows = dest <= remainderNR ? SendNR + 1 : SendNR;
            cout << "Sending " << rows << " rows to task " << dest << " offset= " << offset << endl;
            MPI_Send(&A[offset][0], rows * NCA, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD); // 注意：Send和recv的tag要对应，或者使用同一tag要顺序相同
            MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
            
            MPI_Send(&B, NRB * NCB, MPI_DOUBLE, dest, 4, MPI_COMM_WORLD);
            offset += rows;
        }

        tag = FROM_WORKER;
        for (int dest = 1; dest <= numtasks; dest++)
        { // 接收子进程的计算结果
            MPI_Recv(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, dest, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&C[offset][0], rows * NCB, MPI_DOUBLE, dest, 5, MPI_COMM_WORLD, &status);
            cout << "Received results from task " << dest << endl;
            cout << "worker " << dest << " results is:" << endl;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < NCB; j++)
                {
                    cout << C[offset+i][j] << "\t";
                }
                cout << endl;
            }
        }
        /* 展示计算结果 */
        cout << "C=" << endl;
        for (int i = 0; i < NRA; i++)
        {
            for (int j = 0; j < NCB; j++)
            {
                cout << C[i][j] << "\t";
            }
            cout << endl;
        }
    }

    /* 子进程 */
    if (taskid > MASTER)
    {
        /* 接收主进程的信息 */
        tag = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&A, rows * NCA, MPI_DOUBLE, MASTER, 3, MPI_COMM_WORLD, &status);
        MPI_Recv(&B, NRB * NCB, MPI_DOUBLE, MASTER, 4, MPI_COMM_WORLD, &status);
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < NCB; j++)
            {
                C[i][j] = 0.0;
                for (int k = 0; k < NCA; k++)
                {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        tag = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
        MPI_Send(&C, rows * NCB, MPI_DOUBLE, MASTER, 5, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}