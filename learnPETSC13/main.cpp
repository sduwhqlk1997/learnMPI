#include <petsc.h>
static char help[] = "Mat operator.\n";
int main(int argc, char *argv[])
{
    PetscCall(PetscInitialize(&argc, &argv, NULL, help));
    Mat A, B, C, D, Nest, Dense;
    IS isrow[2], iscol[2];
    Mat mats[4];

    // 创建子矩阵
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
    MatSetFromOptions(A);
    MatSetUp(A);
    MatSetValue(A, 0, 0, 1.0, INSERT_VALUES);
    MatSetValue(A, 1, 1, 1.0, INSERT_VALUES);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    MatCreate(PETSC_COMM_WORLD, &B);
    MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
    //MatSetFromOptions(B);
    MatSetUp(B);
    MatSetValue(B, 0, 1, 2.0, INSERT_VALUES);
    MatSetValue(B, 1, 0, 2.0, INSERT_VALUES);
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

    MatCreate(PETSC_COMM_WORLD, &C);
    MatSetSizes(C, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
    //MatSetFromOptions(C);
    MatSetUp(C);
    MatSetValue(C, 0, 0, 3.0, INSERT_VALUES);
    MatSetValue(C, 1, 1, 3.0, INSERT_VALUES);
    MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);

    MatCreate(PETSC_COMM_WORLD, &D);
    MatSetSizes(D, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
    //MatSetFromOptions(D);
    MatSetUp(D);
    MatSetValue(D, 0, 1, 4.0, INSERT_VALUES);
    MatSetValue(D, 1, 0, 4.0, INSERT_VALUES);
    MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);

    /*PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(MatView(B, PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(MatView(C, PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(MatView(D, PETSC_VIEWER_STDOUT_WORLD));*/

    // 设置嵌套矩阵
    mats[0] = A;
    mats[1] = B;
    mats[2] = C;
    mats[3] = D;

    MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, mats,&Nest);
    MatConvert(Nest,MATDENSE,MAT_INITIAL_MATRIX,&Dense);
    MatView(Dense,PETSC_VIEWER_STDOUT_WORLD);

    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&C);
    MatDestroy(&D);
    MatDestroy(&Nest);

    PetscFinalize();
    return 0;
}