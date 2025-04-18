#include <iostream>
#include <unistd.h>
#include <mpi.h>

#define ARRAYSIZE 20000000
#define MASTER 0

double data[ARRAYSIZE];

double update(int myoffset, int chunk, int myid){
    double mysum;
    mysum = 0;
    for(int i=myoffset; i<myoffset+chunk; i++){
        //data[i] = data[i] + (i*1.0);
        mysum = mysum + data[i];
    }
    printf("Task %d mysum = %e\n",myid,mysum);
    return mysum;
}

int main (int argc, char *argv[]){
   # ifdef DEBUG
    {
        int i=0;
        while (0==i){
            sleep(1);
        }
    }
    #endif

    int numtasks, taskid, rc, dest, offset, tag1, tag2, source, 
        chunksize, leftover;
    double mysum, sum;
    MPI_Status status;
    /*Initializations*/
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    printf ("MPI task %d has started... ",taskid);
    chunksize = ARRAYSIZE/numtasks;
    leftover = ARRAYSIZE % numtasks;
    tag2 = 1;
    tag1 = 2;

    /* Master task only */
    if (taskid == MASTER){
        /* Initialize the array */
        sum = 0;
        for(int i=0;i<ARRAYSIZE;i++){
            data[i] = i*1.0;
            sum+=data[i];
        }
        printf("Initialized array sum = %e\n",sum);
        printf("numtasks = %d chunksize= %d leftover= %d\n",numtasks,chunksize,leftover);
        /* Send each task its portion of the array - master keeps 1st part plus leftover elements */
        offset = chunksize + leftover;
        for (dest=1; dest<numtasks; dest++){
            MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
            MPI_Send(&data[offset], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
            printf("Sent %d elements to task %d offset= %d\n",chunksize,dest,offset);
            offset+=chunksize;
        }
        /* Master does its part of the work */
        offset = 0;
        mysum = update(offset,chunksize+leftover,taskid);

        /* Wait to receive results from each task */
        for (int i=1; i<numtasks;i++){
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            MPI_Recv(&data[offset],chunksize,MPI_DOUBLE,source,tag2,MPI_COMM_WORLD, &status);
        }

        /* Get final sum and print sample results */
        MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

        printf("Sample results: \n");
        offset = 0;
        for (int i=0;i<numtasks;i++){
            for(int j=0;j<5;j++)
                printf(" %e",data[offset+j]);
            printf("\n");
            offset += chunksize;
        }
        printf("*** Final sum= %e ***\n",sum);
    } /* end of master section */

    /* Non-master tasks only */
    if (taskid > MASTER){

        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&data[offset], chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);

        /* Do my part of the work */
        mysum = update(offset, chunksize, taskid);

        /* Send my results back to the master task */
        dest = MASTER;
        MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
        MPI_Send(&data[offset], chunksize, MPI_DOUBLE, MASTER, tag2, MPI_COMM_WORLD);
        /* Use sum reduction operation to obtain final sum */
        MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    } /* end of non-master */

    MPI_Finalize();
}