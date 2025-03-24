#include <iostream>
#include<petsc.h>
#include<cmath>
using namespace std;
int main(int argc, char *argv[]){
    PetscMPIInt rank;
    PetscInt i;
    PetscReal loacalval, globalsum, x;


    PetscCall(PetscInitialize(&argc,&argv,NULL,"compute e in parallel with PETSc.\n\n"));

    PetscOptionsBegin(PETSC_COMM_WORLD,"","options for expx","");
    PetscOptionsReal("-x","input to exp(x) function",NULL,x,&x,NULL);
    PetscOptionsEnd();

    PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
    PetscInt a=1;
    // compute 1/n! where n = (rank of process) + 1
    loacalval = 1.0;
    for(i=2; i<rank+1;i++){
        loacalval/=i;
    }
    for(i=1;i<=rank;i++){
        loacalval*=x;
    }

    //sum the contributions over all process
    PetscCall(MPI_Allreduce(&loacalval, &globalsum,1,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD));

    //output estimate of e and report on work from each process
    PetscPrintf(PETSC_COMM_WORLD,"exp(x) is about %17.15f\n",globalsum);
    PetscReal e,err;
    e=std::exp(x);
    err = abs(globalsum-e)/e;
    PetscPrintf(PETSC_COMM_WORLD,"err=%17.15f\n",err);
    PetscCall(PetscFinalize());
    return 0;
}