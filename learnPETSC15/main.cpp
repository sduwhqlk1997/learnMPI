static char help[] = "A structured-grid Possion solver using DMDA+KSP, FEM2D.\n.\n";
#include "poissonFEM.h"
#include "quadmath.h"
#include <unistd.h>

static PetscReal f_RHS(PetscReal x, PetscReal y){
    return (-y * (1-y)*(1-x-x*x/2) - x*(1-x/2)*(-3*y-y*y)) * PetscExpReal(x+y);
}

static PetscReal uExact(PetscReal x, PetscReal y){
    return x*y*(1-x/2)*(1-y)*PetscExpReal(x+y);
}

int main(int argc, char* argv[]){
    # ifdef DEBUG
    {
        int i=0;
        while (0==i){

            sleep(2);
        }
    }
    #endif
    DM da;
    Mat A;
    Vec b,exact,u;
    KSP ksp;
    PetscReal errnorm, temp;

    DMDALocalInfo info;
    PoissonCtx user;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help));

    user.Lx=1.0;
    user.Ly=1.0;
    user.f_rhs=f_RHS;
    user.g_Dir=uExact;
    user.u_exact=uExact;
    user.quadpts=3;

    // 定义命令
    PetscOptionsBegin(PETSC_COMM_WORLD, "Poisson_","options for main.cpp","");
    PetscCall(PetscOptionsReal("-Lx","","main.cpp", user.Lx,&user.Lx,NULL));
    PetscCall(PetscOptionsReal("-Ly","","main.cpp", user.Ly,&user.Ly,NULL));
    PetscOptionsEnd();

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 
        DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, 3, 3, 
        PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da));
    PetscCall(DMSetApplicationContext(da,&user));
    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));
    PetscCall(DMDASetUniformCoordinates(da,0.0,user.Lx,0.0,user.Ly,0.0,0.0));

    // create Matrix and Vec
    PetscCall(DMCreateGlobalVector(da,&b));
    PetscCall(VecSet(b,0.0));
    PetscCall(DMCreateMatrix(da,&A));
    PetscCall(MatSetFromOptions(A));
    PetscCall(formMatrixVec(da,A,b,&user));
    PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));
    
    // form exact solutions
    PetscCall(DMCreateGlobalVector(da,&exact));
    PetscCall(VecDuplicate(exact,&u));
    PetscCall(formExact(da,exact,&user));

    //create and solve the linear system
    /*
    PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
    PetscCall(KSPSetOperators(ksp,A,A));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp,b,u));
    
    //compute the relative error
    PetscCall(VecAXPY(u,-1.0,exact));
    PetscCall(VecNorm(u,NORM_INFINITY,&errnorm));
    PetscCall(VecNorm(exact,NORM_INFINITY,&temp));
    errnorm /=temp;

    PetscCall(DMDAGetLocalInfo(da,&info));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                          "on %d x %d grid: error |u-uexact|_inf = %g\n",
                        info.mx,info.my,errnorm));
    */
    PetscCall(MatDestroy(&A)); // 销毁矩阵
    //PetscCall(MatDestroy(&B));
    PetscCall(DMDestroy(&da));
    PetscCall(VecDestroy(&b));
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&exact));
    PetscCall(KSPDestroy(&ksp));

    PetscCall(PetscFinalize());
    return 0;
}