#include <iostream>
#include <mpi.h>
int main(int argc,char *argv[]){
    int np, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //定义参数
    int *message=new int[10];
    for(int i=0;i<10;i++){
        message[i]=i+1;
    }

    // 利用reduce函数进行求和
    int *sum=new int[10];
    MPI_Reduce(message,sum,10,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    // 保证所有进程都执行完reduce函数后再执行下面的语句

    if (rank==0){
        for(int i=0;i<10;i++){
            printf("%d ",sum[i]);
        }
        printf("\n");
    }

    // 利用reduce函数进行求积
    int *product = new int[10];
    MPI_Reduce(message,product,10,MPI_INT,MPI_PROD,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0){
        for(int i=0;i<10;i++){
            printf("%d ",product[i]);
        }
        printf("\n");
    }

    // 利用reduce函数求最大值
    int *max = new int[10];
    MPI_Reduce(message,max,10,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0){
        for(int i=0;i<10;i++){
            printf("%d ",max[i]);
        }
        printf("\n");
    }
    MPI_Finalize();

    return 0;
}