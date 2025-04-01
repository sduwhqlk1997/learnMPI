static char help[] = "Solve a tridiagonal system of arbitrary size.\n"
"Option prefix = tri_.\n";

#include <petsc.h>

int main(int argc, char* argv[]){
    Vec x, b, xexact;
    Mat A;
    KSP ksp;
    PetscInt m=4, i, Istart, Iend, j[3];
    PetscReal v[3], xval, errnorm;

    PetscCall(PetscInitialize(&argc,&argv,NULL,help));

    PetscOptionsBegin(PETSC_COMM_WORLD,"tri_","options for tri",NULL);
    PetscCall(PetscOptionsInt("-m","dimension of linear system","main.cpp",m,&m,NULL));
    PetscOptionsEnd();

    PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
    PetscCall(VecSetSizes(x,PETSC_DECIDE,m));
    PetscCall(VecSetFromOptions(x));
    PetscCall(VecDuplicate(x,&b));
    PetscCall(VecDuplicate(x,&xexact));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m));
    PetscCall(MatSetOptionsPrefix(A,"a_")); // 为矩阵 A 的设置命令指定前缀"a_"
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatSetUp(A));
    PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
    for (i=Istart;i<Iend; i++){//设置三对角矩阵
        if (i==0){
            v[0] = 3.0; v[1] = -1.0; // 值
            j[0] = 0; j[1] = 1; // 列
            PetscCall(MatSetValues(A,1,&i,2,j,v,INSERT_VALUES));
        }
        else{
            v[0] = -1.0; v[1] = 3.0; v[2] = -1.0;
            j[0] = i-1; j[1] = i; j[2] = i+1;
            if (i == m - 1){
                PetscCall(MatSetValues(A,1,&i,2,j,v,INSERT_VALUES));
            }
            else{
                PetscCall(MatSetValues(A,1,&i,3,j,v,INSERT_VALUES));
            }
        }
        xval = PetscExpReal(PetscCosReal((double)i)); // e^cos(i)
        PetscCall(VecSetValues(xexact,1,&i,&xval,INSERT_VALUES));
    }
    //从示例代码的第49行开始

}
