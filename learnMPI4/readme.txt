/******************************************************************************
* FILE: mpi_array.c
* DESCRIPTION: 
*   MPI Example - Array Assignment - C Version
*   This program demonstrates a simple data decomposition. The master task
*   first initializes an array and then distributes an equal portion that
*   array to the other tasks. After the other tasks receive their portion
*   of the array, they perform an addition operation to each array element.
*   They also maintain a sum for their portion of the array. The master task 
*   does likewise with its portion of the array and any leftover elements. 
*   As each of the non-master tasks finish, they send their updated portion 
*   of the array to the master.  An MPI collective communication call is used 
*   to collect the local sums maintained by each task.  Finally, the master 
*   task  displays selected parts of the final array and the global sum of 
*   all array elements. 
* AUTHOR: Blaise Barney
* LAST REVISED: 07/03/19
****************************************************************************/