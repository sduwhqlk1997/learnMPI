#ifndef POISSONFEM_H_
#define POISSONFEM_H_
#include <petsc.h>
#include <functional>

namespace{
    using Fun2D = std::function<PetscReal(PetscReal,PetscReal)>;
}

typedef struct 
{
    PetscReal Lx, Ly; // 方程定义区域 x, y 方向的长度
    Fun2D f_rhs;
    Fun2D g_Dir;
    Fun2D u_exact;
}PoissonCtx;

// FEM basis fun
static PetscReal xiL[4] = {1.0, -1.0, -1.0, 1.0},
                 etaL[4] = {1.0, 1.0, -1.0, -1.0}; // 参考单元顶点坐标
typedef struct
{
    PetscInt d_xi, d_eta;
}diff;

static PetscBool Isbound(PetscInt, PetscInt, const DMDALocalInfo*); // 判断当前网格点是否在边界上

static PetscReal chi(PetscInt, PetscReal, diff, PetscReal);

PetscErrorCode formMatrixVec(DM, Mat, Vec, PoissonCtx*);
PetscReal GaussIntegral(Fun2D,PetscReal,PetscReal,PetscInt);
PetscErrorCode formExact(DM, Vec, PoissonCtx*);

#endif