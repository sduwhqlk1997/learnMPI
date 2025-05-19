#include "poissonFEM.h"
#include "quadrature.h"


static PetscReal chi(PetscInt L, PetscReal xi, diff ddiff, PetscReal eta){
    if(ddiff.d_xi==0 && ddiff.d_eta==0)
        return 0.25 * (1.0 + xiL[L] * xi) * (1.0 + etaL[L] * eta);
    else if(ddiff.d_xi==1 && ddiff.d_eta==0){
        return 0.25 * xiL[L] * (1.0 + etaL[L] * eta);
    }
    else if(ddiff.d_xi==0 && ddiff.d_eta==1){
        return 0.25 * etaL[L] * (1.0 + xiL[L] * xi);
    }
    else
        SETERRQ(PETSC_COMM_SELF, 4, "Undefined partial diff");
}

static PetscBool Isbound(PetscInt i, PetscInt j, const DMDALocalInfo* info){
    if(i==0|j==0|i==info->mx-1|j==info->my-1){
        return PETSC_TRUE;
    }
    else{
        return PETSC_FALSE;
    }
}

PetscErrorCode formMatrixVec(DM da, Mat A, Vec b, PoissonCtx *user){
    DMDALocalInfo info;
    MatStencil Arow, Acol;
    const PetscInt li[4] = {0,1,1,0}, lj[4] = {0,0,1,1};
    for(PetscInt j=info.ys; j<info.ys+info.ym;j++){
        if(j==info.my-1)
            continue;
        for(PetscInt i=info.xs; i<info.xs+info.xm;i++){
            if(i==info.mx-1)
                continue;
            for(PetscInt l=0; l<4; l++){//测试函数(行)
                PetscInt itest = i+li[l], jtest=j+lj[l];
                if(Isbound(itest,jtest,&info)==PETSC_TRUE){//在边界上
                    for(PetscInt r=0; r<4; r++){//试探函数(列)

                    }
                }else{
                    for(PetscInt r=0; r<4; r++){//试探函数(列)
                        PetscInt itrail = i+li[r], jtrail=j+lj[r];
                        if(Isbound(itrail,jtrail,&info)==PETSC_TRUE){//在边界上

                        }else{//test和trail全为内部点

                        }
                    }
                }
            }
        }
    }
    return PETSC_SUCCESS;
}

PetscReal GaussIntegral(Fun2D IntFun, PetscReal hx, PetscReal hy, PetscInt order)
{
    const Quad1D q = gausslegendre[order - 1];
    PetscReal Int = 0;
    PetscReal Jacobi = 0.25 * hx * hy;
    for (PetscInt i = 0; i < order; i++)
    {
        for (PetscInt j = 0; j < order; j++)
        {
            Int += q.w[i] * q.w[j] * IntFun(q.xi[i], q.xi[j]);
        }
    }
    return Jacobi * Int;
}

PetscErrorCode formExact(DM da, Vec uexact, PoissonCtx * user){
    PetscReal hx, hy, x, y, **auexact, xymin[2], xymax[2];
    DMDALocalInfo info;
    PetscCall(DMDAGetLocalInfo(da,&info));
    PetscCall(DMGetBoundingBox(info.da,xymin,xymax));

    hx = user->Lx/(info.mx-1); hy = user->Ly/(info.my-1);
    PetscCall(DMDAVecGetArray(da, uexact, &auexact));
    for(PetscInt j = info.ys; j<info.ys + info.ym; j++){
        y = xymin[1]+hy*j;
        for(PetscInt i = info.xs;i<info.xs+info.xm;i++){
            x = xymin[0]+hx*i;
            auexact[j][i] = user->u_exact(x,y);
        }
    }
    PetscCall(DMDAVecRestoreArray(da,uexact,&auexact));
    return PETSC_SUCCESS;
}