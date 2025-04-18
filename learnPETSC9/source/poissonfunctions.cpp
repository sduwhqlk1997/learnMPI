# include <petsc.h>
#include "poissonfunctions.h"

PetscErrorCode Poisson1DFunctionLocal(DMDALocalInfo *info, PetscReal *au,
    PetscReal *aF, PoissonCtx *user) {
PetscInt   i;
PetscReal  xmax[1], xmin[1], h, x, ue, uw;
PetscCall(DMGetBoundingBox(info->da,xmin,xmax)); // 获取网格的边界
h = (xmax[0] - xmin[0]) / (info->mx - 1);
for (i = info->xs; i < info->xs + info->xm; i++) {
x = xmin[0] + i * h;
if (i==0 || i==info->mx-1) {
aF[i] = au[i] - user->g_bdry(x,0.0,0.0,user);
aF[i] *= user->cx * (2.0 / h);
} else {
ue = (i+1 == info->mx-1) ? user->g_bdry(x+h,0.0,0.0,user)
   : au[i+1];
uw = (i-1 == 0)          ? user->g_bdry(x-h,0.0,0.0,user)
   : au[i-1];
aF[i] = user->cx * (2.0 * au[i] - uw - ue) / h
- h * user->f_rhs(x,0.0,0.0,user);
}
}
PetscCall(PetscLogFlops(9.0*info->xm)); // 计算浮点运算个数
return 0;
}

//STARTFROM2DFUNCTION
PetscErrorCode Poisson2DFunctionLocal(DMDALocalInfo *info, PetscReal **au, PetscReal **aF, PoissonCtx *user){
    PetscInt i,j;
    PetscReal xymin[2], xymax[2], hx, hy, darea, scx, scy, scdiag, x, y, ue, uw, un, us;
    PetscCall(DMGetBoundingBox(info->da,xymin,xymax));
    hx = (xymax[0]-xymin[0]) / (info->mx - 1);
    hy = (xymax[1]-xymin[1]) / (info->my - 1);
    darea = hx * hy;
    scx = user->cx * hy / hx;
    scy = user->cy * hx / hy;
    scdiag = 2.0 * (scx + scy); // diagonal scaling
    for (j = info->ys; j<info->ys + info->ym; j++){
        y = xymin[1] + j * hy;
        for(i = info->xs; i<info->xs + info->xm; i++){
            x = xymin[0] + i * hx;
            if (i==0 || i==info->mx-1 || j==0 || j==info->my-1){
                aF[j][i] = au[j][i] - user->g_bdry(x,y,0.0,user);
                aF[j][i] *= scdiag;
            } else {
                ue = (i+1 == info->mx-1) ? user->g_bdry(x+hx,y,0.0,user) : au[j][i+1];
                uw = (i-1 == 0) ? user->g_bdry(x-hx,y,0.0,user) : au[j][i-1];
                un = (j+1 == info->my-1) ? user->g_bdry(x,y+hy,0.0,user) : au[j+1][i];
                us = (j-1 == 0) ? user->g_bdry(x,y-hy,0.0,user) : au[j-1][i];
                aF[j][i] = scdiag * au[j][i] - scx * (uw + ue) - scy * (us + un) 
                            - darea * user->f_rhs(x,y,0.0,user);
            }
        }
    }
    PetscCall(PetscLogFlops(11.0*info->xm*info->ym));
    return 0;
}
//ENDFORM2DFUNCTION

PetscErrorCode Poisson3DFunctionLocal(DMDALocalInfo *info, PetscReal ***au,
    PetscReal ***aF, PoissonCtx *user) {
PetscInt   i, j, k;
PetscReal  xyzmin[3], xyzmax[3], hx, hy, hz, dvol, scx, scy, scz, scdiag,
x, y, z, ue, uw, un, us, uu, ud;
PetscCall(DMGetBoundingBox(info->da,xyzmin,xyzmax));
hx = (xyzmax[0] - xyzmin[0]) / (info->mx - 1);
hy = (xyzmax[1] - xyzmin[1]) / (info->my - 1);
hz = (xyzmax[2] - xyzmin[2]) / (info->mz - 1);
dvol = hx * hy * hz;
scx = user->cx * dvol / (hx*hx);
scy = user->cy * dvol / (hy*hy);
scz = user->cz * dvol / (hz*hz);
scdiag = 2.0 * (scx + scy + scz);
for (k = info->zs; k < info->zs + info->zm; k++) {
z = xyzmin[2] + k * hz;
for (j = info->ys; j < info->ys + info->ym; j++) {
y = xyzmin[1] + j * hy;
for (i = info->xs; i < info->xs + info->xm; i++) {
x = xyzmin[0] + i * hx;
if (   i==0 || i==info->mx-1
|| j==0 || j==info->my-1
|| k==0 || k==info->mz-1) {
aF[k][j][i] = au[k][j][i] - user->g_bdry(x,y,z,user);
aF[k][j][i] *= scdiag;
} else {
ue = (i+1 == info->mx-1) ? user->g_bdry(x+hx,y,z,user)
           : au[k][j][i+1];
uw = (i-1 == 0)          ? user->g_bdry(x-hx,y,z,user)
           : au[k][j][i-1];
un = (j+1 == info->my-1) ? user->g_bdry(x,y+hy,z,user)
           : au[k][j+1][i];
us = (j-1 == 0)          ? user->g_bdry(x,y-hy,z,user)
           : au[k][j-1][i];
uu = (k+1 == info->mz-1) ? user->g_bdry(x,y,z+hz,user)
           : au[k+1][j][i];
ud = (k-1 == 0)          ? user->g_bdry(x,y,z-hz,user)
           : au[k-1][j][i];
aF[k][j][i] = scdiag * au[k][j][i]
- scx * (uw + ue) - scy * (us + un) - scz * (uu + ud)
- dvol * user->f_rhs(x,y,z,user);
}
}
}
}
PetscCall(PetscLogFlops(14.0*info->xm*info->ym*info->zm));
return 0;
}