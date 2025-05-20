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
    const PetscInt li[4] = {0,1,1,0}, lj[4] = {0,0,1,1};
    PetscReal hx,hy,x,y,kx,ky,tx,ty,cx,cy, xymin[2], xymax[2], xEleStart, yEleStart;
    PetscCall(DMDAGetLocalInfo(da,&info));
    PetscCall(DMGetBoundingBox(info.da,xymin,xymax));

    hx = user->Lx/(info.mx-1); hy = user->Ly/(info.my-1);
    ky = hy / 2.0; kx = hx / 2.0;
    cx = 4/(hx*hx); cy = 4/(hy*hy);

    MatStencil* AIndexInsert = new MatStencil[info.mx*info.my];
    PetscReal* AValueInsert = new PetscReal[info.mx*info.my];

    MatStencil* ArowAdd = new MatStencil[4*info.mx*info.my];
    MatStencil* AcolAdd = new MatStencil[4*info.mx*info.my];
    PetscReal* AValueAdd = new PetscReal[4*info.mx*info.my];

    PetscInt* bindexInsert = new PetscInt[info.mx*info.my];
    PetscReal* bValueInsert = new PetscReal[info.mx*info.my];
    
    PetscInt* bindexAdd = new PetscInt[4*info.mx*info.my];
    PetscReal* bValueAdd = new PetscReal[4*info.mx*info.my];

    PetscInt ANInsert = 0, ANAdd = 0, bNInsert=0, bNAdd=0;

    for(PetscInt j=info.ys; j<info.ys+info.ym;j++){
        if(j==info.my-1)
            continue;
        for(PetscInt i=info.xs; i<info.xs+info.xm;i++){
            if(i==info.mx-1)
                continue;
            tx = xymin[0] + i * hx + kx; //仿射变换的bias
            ty = xymin[1] + j * hy + ky;
            for(PetscInt l=0; l<4; l++){//测试函数(行)
                PetscInt itest = i+li[l], jtest=j+lj[l];
                if(Isbound(itest,jtest,&info)==PETSC_TRUE){//边界点对应的行
                    // 矩阵
                    AIndexInsert[ANInsert].i=itest; AIndexInsert[ANInsert].j=jtest;
                    AIndexInsert[ANInsert].k=0; AIndexInsert[ANInsert].c=0;
                    AValueAdd[ANInsert++] = 1.0;
                    // 向量 直接将边界上的解赋值
                    x = xymin[0]+hx*itest; y = xymin[1]+hy*jtest; // 网格点坐标
                    bindexInsert[bNInsert] = info.mx*jtest+itest;
                    bValueInsert[bNInsert++] = user->g_Dir(x,y);
                }else{
                    // 向量 此时测试函数为内部点，做积分
                    Fun2D IntFunVec=[l,kx,ky,tx,ty,user](PetscReal x, PetscReal y){
                        return chi(l,x,{0,0},y) * user->f_rhs(kx*x+tx,ky*y+ty);
                    };
                    bindexAdd[bNAdd] = info.mx*jtest+itest;
                    bValueAdd[bNAdd++] = GaussIntegral(IntFunVec, hx, hy, user->quadpts);

                    for(PetscInt r=0; r<4; r++){//试探函数(列)
                        PetscInt itrail = i+li[r], jtrail=j+lj[r];
                        Fun2D IntFunMat = [l,r,cx,cy](PetscReal x,PetscReal y){
                                return chi(l,x,{1,0},y)*chi(r,x,{1,0},y)*cx+
                                chi(l,x,{0,1},y)*chi(r,x,{0,1},y)*cy;
                        };
                        if(Isbound(itrail,jtrail,&info)==PETSC_TRUE){//边界点对应的列
                            x = xymin[0] + itrail*hx; y = xymin[1] + jtrail*hy;
                            PetscReal boundVal = user->g_Dir(x,y);
                            bindexAdd[bNAdd] = info.mx*jtest + itest;
                            bValueAdd[bNAdd++] = -boundVal*GaussIntegral(IntFunMat,hx,hy,user->quadpts);
                        }else{//test和trail全为内部点
                            ArowAdd[ANAdd].c=0; ArowAdd[ANAdd].k=0; ArowAdd[ANAdd].i=itest; ArowAdd[ANAdd].j=jtest;
                            AcolAdd[ANAdd].c=0; AcolAdd[ANAdd].k=0; AcolAdd[ANAdd].i=itrail; AcolAdd[ANAdd].j=jtrail;
                            AValueAdd[ANAdd++] = GaussIntegral(IntFunMat,hx,hy,user->quadpts);
                        }
                    }
                }
            }
        }
    }
    // 向矩阵和向量中填充值
    for(PetscInt i=0;i<ANAdd;i++){
        PetscCall(MatSetValuesStencil(A,1,&ArowAdd[i],1,&AcolAdd[i],&AValueAdd[i],ADD_VALUES));
    }
    PetscCall(MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY));

    for(PetscInt i=0;i<ANInsert;i++){
        PetscCall(MatSetValuesStencil(A,1,&AIndexInsert[i],1,&AIndexInsert[i],&AValueInsert[i],INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

    PetscCall(VecSetValues(b,bNAdd,bindexAdd,bValueAdd,ADD_VALUES));
    PetscCall(VecAssemblyBegin(b));
    PetscCall(VecAssemblyEnd(b));

    PetscCall(VecSetValues(b,bNInsert,bindexInsert,bValueInsert,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(b));
    PetscCall(VecAssemblyEnd(b));
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