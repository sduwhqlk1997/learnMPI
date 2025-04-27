static char help[] = "A structured-grid Possion solver using DMDA+KSP, FEM2D.\n.\n";
#include <petsc.h>
#include <functional>
#include <unistd.h>
#include "quadrature.h"

typedef struct
{
    PetscInt quadpts;
}PHelmCtx;

using integrandFun = std::function<PetscReal(PetscReal,PetscReal)>;
//typedef PetscReal (*integrandFun)(PetscReal,PetscReal); // 被积函数
//PetscReal GaussIntegral(std::function<PetscReal(PetscReal,PetscReal)>,PetscReal*,PetscReal*,PetscInt); // 数值积分
PetscReal GaussIntegral(integrandFun,
    PetscReal,PetscReal,
    PetscInt); // 数值积分
PetscErrorCode formMatrix(DM, Mat,PHelmCtx*);
PetscErrorCode formRHS(DM, Vec);
PetscErrorCode formExact(DM, Vec);

// FEM basis fun
static PetscReal xiL[4] = {1.0, -1.0, -1.0, 1.0},
                 etaL[4] = {1.0, 1.0, -1.0, -1.0}; // 参考单元顶点坐标
typedef struct
{
    PetscInt d_xi, d_eta;
}diff;

static PetscReal chi(PetscInt, PetscReal, diff, PetscReal);

int main(int argc, char* argv[]){
    # ifdef DEBUG
    {
        int i=0;
        while (0==i){
            sleep(10);
        }
    }
    #endif
    DM da;
    Mat A;
    DMDALocalInfo info;
    PHelmCtx user;
    user.quadpts = 3;

    PetscCall(PetscInitialize(&argc,&argv,NULL,help));

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 
        DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, 3, 3, 
        PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da));

    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));
    PetscCall(DMCreateMatrix(da,&A));
    //MatSetType(A, MATSBAIJ);
    PetscCall(MatSetFromOptions(A));
    PetscCall(formMatrix(da,A,&user));
    PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));

    PetscCall(DMDAGetLocalInfo(da,&info));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "done on %d x %d grid...\n",
        info.mx, info.my));
    
    PetscCall(MatDestroy(&A)); // 销毁矩阵
    PetscCall(DMDestroy(&da));
    PetscCall(PetscFinalize());
    return 0;
}

PetscReal GaussIntegral(integrandFun IntFun,
    PetscReal hx,PetscReal hy,
    PetscInt order){
        const Quad1D q = gausslegendre[order-1];
        PetscReal Int = 0;
        PetscReal Jacobi = 0.25*hx*hy;
        for(PetscInt i = 0; i<order; i++){
            for(PetscInt j = 0; j<order; j++){
                Int += q.w[i] * q.w[j] * IntFun(q.xi[i],q.xi[j]);
            }
        }
        return Jacobi * Int;
    }

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

PetscErrorCode formMatrix(DM da, Mat A, PHelmCtx* user){ // A 为对称矩阵，MatSetType(A, MATSBAIJ)
    DMDALocalInfo info;
    MatStencil row, col;
    PetscReal hx, hy, v;
    const PetscInt li[4] = {0,1,1,0}, lj[4] = {0,0,1,1};
    PetscInt nrows, ncols;

    PetscCall(DMDAGetLocalInfo(da,&info));
    hx = 1.0/(info.mx-1); hy = 1.0/(info.my-1);
    for(PetscInt j=info.ys; j<info.ys+info.ym; j++){// 遍历定点为(i,j),(i+1,j),(i+1,j+1),(i,j+1)的单元
        if(j==info.my-1)
            continue;
        for(PetscInt i=info.xs; i<info.xs+info.xm; i++){
            if(i==info.mx-1)
                continue;
            PetscReal cx = 4/(hx*hx), cy = 4/(hy*hy);
            for(PetscInt l=0; l<4 ;l++){// 遍历每个单元的基函数
                for(PetscInt r = 0; r<4; r++){
                    //PetscInt index=0;
                    integrandFun IntFun=[l,r,cx,cy](PetscReal x,PetscReal y){
                        return chi(l,x,{1,0},y)*chi(r,x,{1,0},y)*cx+
                        chi(l,x,{0,1},y)*chi(r,x,{0,1},y)*cy;
                    };
                    v=GaussIntegral(IntFun,hx,hy,user->quadpts);
                    row.i=i+li[l]; row.j=j+lj[l]; row.k = 0; row.c = 0;// 测试函数
                    col.i=i+li[r]; col.j=j+lj[r]; col.k - 0; col.c = 0;// 试探函数
                    PetscCall(MatSetValuesStencil(A,1,&row,1,&col,&v,ADD_VALUES));
                }
            }
        }
    }
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
    return PETSC_SUCCESS;
}

