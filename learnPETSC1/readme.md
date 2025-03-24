PETSC for PDE 的算例1
逼近自然数e
#########笔记##########
1. 在include<petsc.h>后实际上已经包含了mpi.h
2. PetscReal数据类型与double类型是相似的
3. 使用PetscPrintf()打印时，若 communicator 设置为 <PETSC_COMM_WORLD> 则只会在主线程(rank zero)打印一次；若使用 <PETSC_COMM_SELF> 则会在每个线程都打印一次；若要进行有序输出，则可调用<PetscSynchronizedPrintf()>，不过这需要额外的进程间通信。
