#include <iostream>
#include<petsc.h>
#include<cmath>
using namespace std;
int main(int argc, char *argv[]){
    PetscMPIInt rank;
    PetscInt i;
    PetscReal loacalval, globalsum;

    PetscCall(PetscInitialize(&argc,&argv,NULL,"compute e in parallel with PETSc.\n\n"));
    PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
    PetscInt a=1;
    // compute 1/n! where n = (rank of process) + 1
    loacalval = 1.0;
    for(i=2; i<rank+1;i++){
        loacalval/=i;
    }

    //sum the contributions over all process
    PetscCall(MPI_Allreduce(&loacalval, &globalsum,1,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD));

    //output estimate of e and report on work from each process
    PetscPrintf(PETSC_COMM_WORLD,"e is about %17.15f\n",globalsum);
    PetscReal e,err;
    e=2.71828182845904523536;
    err = abs(globalsum-e)/e;
    PetscPrintf(PETSC_COMM_WORLD,"err=%17.15f\n",err);
    PetscCall(PetscFinalize());
    MPI_Finalize();
    cout<<"a"<<endl;
    return 0;
}