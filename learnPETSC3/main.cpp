#include <petsc.h>
static char help[] = "Assemble a Mat sparsely.\n"; // 若不加关键字 "static"，在其他文件也有 help变量时会有冲突
int main(int argc,char *argv[]){
    Mat A;
    PetscInt i1[3] = {0,1,2},
             j1[3] = {0,1,2},
             i2 = 3,
             j2[3] = {1,2,3},
             i3 = 1,
             j3 = 3;
    PetscReal aA1[9] = {1.0, 2.0, 3.0,
                        2.0, 1.0, -2.0,
                        -1.0, 1.0, 1.0},
              aA2[3] = {1.0, 1.0, -1.0},
              aA3 = -3.0;
    
    PetscCall(PetscInitialize(&argc,&argv,NULL,help));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,4,4));
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatSetUp(A)); // 初始化矩阵

    PetscCall(MatSetValues(A,3,i1,3,j1,aA1,INSERT_VALUES));
    PetscCall(MatSetValues(A,1,&i2,3,j2,aA2,INSERT_VALUES)); // 填入多个值用MatSetValues
    PetscCall(MatSetValue(A,i3,j3,aA3,INSERT_VALUES)); // 只填一个值用MatSetValue

    /*完成矩阵的组装过程，确保矩阵的值和结构已经准备好用于后续计算
    * MAT_FINAL_ASSEMBLY：表示这是最终的组装阶段（即矩阵的值和结构将不再更改）
    * 如果是并行矩阵，MatAssemblyBegin 会开始进程间的通信，以收集其他进程对矩阵的修改。
    * 如果是并行矩阵，MatAssemblyEnd 会结束进程间的通信，并确保所有进程的修改都已同步。
    */
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY)); // 完成矩阵的组装过程，确保矩阵的值和结构已经准备好用于后续计算

    PetscCall(MatDestroy(&A)); // 销毁矩阵
    PetscCall(PetscFinalize());
    return 0;
}