static char help[] = "delete some rows and cols.\n";
#include <petsc.h>

int main(int argc, char* argv[]){
    PetscInitialize(&argc, &argv, NULL, help);
    Mat A; // 原矩阵
    Mat subA; // 子矩阵
    IS isRow, isCol; // 定义要保留的行和列的索引集合
    PetscInt rowsToKeep[] = {0, 2, 3}; // 保留的行索引
    PetscInt colsToKeep[] = {1, 3}; // 保留的列索引
    PetscInt nRows = 3, nCols = 2;

    // 创建矩阵 A（4x4 矩阵）
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 4, 4); // 4x4 矩阵
    MatSetFromOptions(A);
    MatSetUp(A);

    // 按行填充矩阵 A 为 1, 2, 3, ..., 16
    PetscInt i, j;
    PetscScalar value = 1.0;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            MatSetValue(A, i, j, value, INSERT_VALUES);
            value += 1.0;
        }
    }

    // 完成矩阵组装
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // 查看原矩阵 A
    PetscPrintf(PETSC_COMM_WORLD, "Original Matrix A:\n");
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // 创建行和列的索引集合
    ISCreateGeneral(PETSC_COMM_WORLD, nRows, rowsToKeep, PETSC_COPY_VALUES, &isRow);
    ISCreateGeneral(PETSC_COMM_WORLD, nCols, colsToKeep, PETSC_COPY_VALUES, &isCol);

    // 创建子矩阵，排除不需要的行和列
    MatCreateSubMatrix(A, isRow, isCol, MAT_INITIAL_MATRIX, &subA);

    PetscPrintf(PETSC_COMM_WORLD, "\nSubmatrix after removing rows and columns:\n");
    MatView(subA, PETSC_VIEWER_STDOUT_WORLD);

    // 销毁对象
    ISDestroy(&isRow);
    ISDestroy(&isCol);
    MatDestroy(&A);
    MatDestroy(&subA);

    PetscFinalize();
    return 0;
}