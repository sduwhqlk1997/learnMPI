static char help[] = "Unstructured 2D FEM solution of nolinear Poisson problem. Option prefix un_.\n";

#include "quadrature.h"
#include "um.h"
#include "case.h"

// STARTCTX
typedef struct{
    UM *mesh;
    PetscInt solncase,
             quaddegree;
    PetscReal (*a_fcn)(PetscReal,PetscReal,PetscReal);
    PetscReal (*f_fcn)(PetscReal,PetscReal,PetscReal);
    PetscReal (*gD_fcn)(PetscReal,PetscReal);
    PetscReal (*gN_fcn)(PetscReal,PetscReal);
    PetscReal (*uexact_fcn)(PetscReal,PetscReal);
    PetscLogStage readstage, setupstage,solverstage,resstage,jacstage;
}umfemCtx;
// ENDCTX

//STARTFEM
PetscReal chi(PetscInt L,PetscReal xi, PetscReal eta){
    const PetscReal z[3] = {1.0 - xi - eta, xi, eta};
    return z[L];
}

const PetscReal dchi[3][2] = {{-1.0,-1.0},{1.0,0.0},{0.0,1.0}};

// evaluate v(xi,eta) on reference element using local node numbering
PetscReal eval(const PetscReal v[3],PetscReal xi,PetscReal eta){
    PetscReal sum = 0.0;
    PetscInt L;
    for(L=0;L<3;L++)
        sum+=v[L]*chi(L,xi,eta);
    return sum;
}
//ENDFEM

extern PetscErrorCode FillExact(Vec,umfemCtx*);
extern PetscErrorCode FormFunction(SNES, Vec, Vec, void*);
extern PetscErrorCode FormPicard(SNES, Vec, Mat, Mat, void*);
extern PetscErrorCode preallocateAndSetNonzeros(Mat, umfemCtx*);

int main(int argc, char* argv[]){
    PetscMPIInt size;
    PetscBool viewmesh = PETSC_FALSE,
              viewsoln = PETSC_FALSE,
              noprealloc = PETSC_FALSE,
              savepintbinary = PETSC_FALSE,
              savepintmatlab = PETSC_FALSE;
    char root[256] = "", nodesname[256], issname[256], solnname[256], pintname[256] = "";
    PetscInt savepintlevel = -1, levels;
    UM mesh;
    umfemCtx user;
    SNES snes;
    KSP ksp;
    PC pc;
    PCType pctype;
    Mat A;
    Vec r,u,uexact;
    PetscReal err, h_max;

    PetscCall(PetscInitialize(&argc,&argv,NULL,help));

    PetscCall(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    if(size!=1){
        SETERRQ(PETSC_COMM_SELF,1,"umfem only works on one MPI process");
    }


}