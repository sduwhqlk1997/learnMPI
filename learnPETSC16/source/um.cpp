#include "um.h"

PetscErrorCode UMInitialize(UM *mesh){
    mesh->N = 0;
    mesh->K = 0;
    mesh->P = 0;
    mesh->loc = NULL;
    mesh->e = NULL;
    mesh->bf = NULL;
    mesh->ns = NULL;
    return 0;
}

PetscErrorCode UMDestory(UM *mesh){
    PetscCall(VecDestroy(&(mesh->loc)));
    PetscCall(ISDestroy(&(mesh->e)));
    PetscCall(ISDestroy(&(mesh->bf)));
    PetscCall(ISDestroy(&(mesh->ns)));
    return 0;
}

PetscErrorCode UMReadNodes(UM *mesh, char *filename){
    PetscInt twoN;
    PetscViewer viewer;
    if (mesh->N > 0){
        SETERRQ(PETSC_COMM_SELF,1,"nodes already created?\n");
    }
    PetscCall(VecCreate(PETSC_COMM_WORLD,&mesh->loc));
    PetscCall(VecSetFromOptions(mesh->loc));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
    PetscCall(VecLoad(mesh->loc,viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(VecGetSize(mesh->loc,&twoN));
    if (twoN % 2 != 0){
        SETERRQ(PETSC_COMM_SELF,2,"node locations loaded from %s are not N pairs\n",filename);
    }
    mesh->N = twoN / 2;
    return 0;
}

PetscErrorCode UMReadISs(UM *mesh, char *filename){
    PetscViewer viewer;
    PetscInt n_bf;
    if ((!mesh->loc) || (mesh->N==0)){
        SETERRQ(PETSC_COMM_SELF,2,"node cooedinates not created ... do that first ... stopping\n");
    }
    if ((mesh->K>0) || (mesh->P>0) || (mesh->e!=NULL) || (mesh->bf!=NULL)||(mesh->ns != NULL)){
        SETERRQ(PETSC_COMM_SELF,1,"elements, boundary flags, Neumann boundary segments already created?... stopping\n");
    }
    return 0;
}