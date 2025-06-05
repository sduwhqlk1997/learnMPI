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
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
    // create and load e
    PetscCall(ISCreate(PETSC_COMM_WORLD,&(mesh->e)));
    PetscCall(ISLoad(mesh->e,viewer));
    PetscCall(ISGetSize(mesh->e,&(mesh->K)));
    if(mesh->K % 3 != 0){
        SETERRQ(PETSC_COMM_SELF,3,"IS e loaded from %s is wrong size for list of element triples\n",filename);
    }
    mesh->K /= 3;
    // create and load bf
    PetscCall(ISCreate(PETSC_COMM_WORLD,&(mesh->bf)));
    PetscCall(ISLoad(mesh->bf,viewer));
    PetscCall(ISGetSize(mesh->bf,&n_bf));
    if(n_bf != mesh->N){
        SETERRQ(PETSC_COMM_SELF,4,"IS bf loaded from %s is wrong size for list of boundary flags\n",filename);
    }
    // FIXME seems there is no way to tell if file is empty at this point
    // create and load ns last ... may *start with a negative value* in which case set P=0
    const PetscInt *ans;
    PetscCall(ISCreate(PETSC_COMM_WORLD,&(mesh->ns)));
    PetscCall(ISLoad(mesh->ns,viewer));
    PetscCall(ISGetIndices(mesh->ns,&ans));
    if (ans[0]<0){
        ISDestroy(&(mesh->ns));
        mesh->ns = NULL;
        mesh->P = 0;
    } else{
        PetscCall(ISGetSize(mesh->ns,&(mesh->P)));
        if(mesh->P % 2 != 0){
            SETERRQ(PETSC_COMM_SELF,4,"IS s loaded from %s is wrong size for list of Neumann boundary segment pairs\n",filename);
        }
        mesh->P /= 2;
    }
    PetscCall(PetscViewerDestroy(&viewer));

    // check that mesh is complete now
    PetscCall(UMCheckElements(mesh));
    PetscCall(UMCheckBoundaryData(mesh));
    return 0;
}

PetscErrorCode UMCheckElements(UM *mesh){
    const PetscInt *ae;
    PetscInt k,m;
    if ((mesh->K == 0) || (mesh->e == NULL)){
        SETERRQ(PETSC_COMM_SELF,1,"number of elements unknown; call UMReadElements() first\n");
    }
    if(mesh->N == 0){
        SETERRQ(PETSC_COMM_SELF,2,"node size unknown so element check impossible; call UMReadNodes() first\n");
    }
    PetscCall(ISGetIndices(mesh->e,&ae));
    for (k=0;k<mesh->K;k++){
        for(m=0;m<3;m++){
            if ((ae[3*k+m]<0) || (ae[3*k+m]>=mesh->N)){
                SETERRQ(PETSC_COMM_SELF,3,"index e[%d]=%d invalid: not between 0 and N-1=%d\n",3*k+m,ae[3*k+m],mesh->N-1);
            }
        }
    }
    PetscCall(ISRestoreIndices(mesh->e,&ae));
    return 0;
}

PetscErrorCode UMCheckBoundaryData(UM *mesh){
    const PetscInt *ans, *abf;
    PetscInt n,m;
    if(mesh->N == 0){
        SETERRQ(PETSC_COMM_SELF,2,"node size unknown so boundary flag check impossible; call UMReadNodes() first\n");
    }
    if(mesh->bf == NULL){
        SETERRQ(PETSC_COMM_SELF,1,"boundary flags at nodes not allocated; call UMReadNodes() first\n");
    }
    if((mesh->P>0)&&(mesh->ns==NULL)){
        SETERRQ(PETSC_COMM_SELF,3,"inconsistent data for Neumann boundary segments\n");
    }
    PetscCall(ISGetIndices(mesh->bf,&abf));
    for(n=0;n<mesh->N; n++){
        switch(abf[n]){
            case 0:
            case 1:
            case 2:
                break;
            default :
                SETERRQ(PETSC_COMM_SELF,5,"boundary flag bf[%d]=%d invalid: not in {0,1,2}\n",n,abf[n]);
        }
    }
    PetscCall(ISRestoreIndices(mesh->bf,&abf));
    if(mesh->P > 0){
        PetscCall(ISGetIndices(mesh->ns,&ans));
        for(n=0;n<mesh->P;n++){
            for (m=0;m<2;m++){
                if ((ans[2*n+m]<0) || (ans[2*n+m]>=mesh->N)){
                    SETERRQ(PETSC_COMM_SELF,6,"index ns[%d]=%d invalid: not between 0 and N-1=%d\n",2*n+m,ans[3*n+m],mesh->N-1);
                }
            }
        }
        PetscCall(ISRestoreIndices(mesh->ns,&ans));
    }
    return 0;
}