#include <petsc.h>
#include <iostream>

int main(int argc, char **argv) {
    // 初始化 PETSc
    PetscInitialize(&argc, &argv, (char *)0, NULL);

    // 获取 PETSc 版本信息
    char *version = new char[1000];
    size_t int_size;
    PetscGetVersion(version,int_size); // 传入 version 的地址

    // 输出 PETSc 版本信息
    std::cout << "PETSc version: " << version << std::endl;

    // 获取处理器数量
    int size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    // 输出处理器数量
    std::cout << "Number of processors: " << size << std::endl;

    // 结束 PETSc
    PetscFinalize();

    return 0;
}