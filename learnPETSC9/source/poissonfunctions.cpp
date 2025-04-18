#include <petsc.h>
#include "poissonfunctions.h"

PetscErrorCode Poisson1DFunctionLocal(DMDALocalInfo *info, PetscReal *au,
                                      PetscReal *aF, PoissonCtx *user)
{
    PetscInt i;
    PetscReal xmax[1], xmin[1], h, x, ue, uw;
    PetscCall(DMGetBoundingBox(info->da, xmin, xmax)); // 获取网格的边界
    h = (xmax[0] - xmin[0]) / (info->mx - 1);
    for (i = info->xs; i < info->xs + info->xm; i++)
    {
        x = xmin[0] + i * h;
        if (i == 0 || i == info->mx - 1)
        {
            aF[i] = au[i] - user->g_bdry(x, 0.0, 0.0, user);
            aF[i] *= user->cx * (2.0 / h);
        }
        else
        {
            ue = (i + 1 == info->mx - 1) ? user->g_bdry(x + h, 0.0, 0.0, user)
                                         : au[i + 1];
            uw = (i - 1 == 0) ? user->g_bdry(x - h, 0.0, 0.0, user)
                              : au[i - 1];
            aF[i] = user->cx * (2.0 * au[i] - uw - ue) / h - h * user->f_rhs(x, 0.0, 0.0, user);
        }
    }
    PetscCall(PetscLogFlops(9.0 * info->xm)); // 计算浮点运算个数
    return 0;
}

// STARTFROM2DFUNCTION
PetscErrorCode Poisson2DFunctionLocal(DMDALocalInfo *info, PetscReal **au, PetscReal **aF, PoissonCtx *user)
{
    PetscInt i, j;
    PetscReal xymin[2], xymax[2], hx, hy, darea, scx, scy, scdiag, x, y, ue, uw, un, us;
    PetscCall(DMGetBoundingBox(info->da, xymin, xymax));
    hx = (xymax[0] - xymin[0]) / (info->mx - 1);
    hy = (xymax[1] - xymin[1]) / (info->my - 1);
    darea = hx * hy;
    scx = user->cx * hy / hx;
    scy = user->cy * hx / hy;
    scdiag = 2.0 * (scx + scy); // diagonal scaling
    for (j = info->ys; j < info->ys + info->ym; j++)
    {
        y = xymin[1] + j * hy;
        for (i = info->xs; i < info->xs + info->xm; i++)
        {
            x = xymin[0] + i * hx;
            if (i == 0 || i == info->mx - 1 || j == 0 || j == info->my - 1)
            {
                aF[j][i] = au[j][i] - user->g_bdry(x, y, 0.0, user);
                aF[j][i] *= scdiag;
            }
            else
            {
                ue = (i + 1 == info->mx - 1) ? user->g_bdry(x + hx, y, 0.0, user) : au[j][i + 1];
                uw = (i - 1 == 0) ? user->g_bdry(x - hx, y, 0.0, user) : au[j][i - 1];
                un = (j + 1 == info->my - 1) ? user->g_bdry(x, y + hy, 0.0, user) : au[j + 1][i];
                us = (j - 1 == 0) ? user->g_bdry(x, y - hy, 0.0, user) : au[j - 1][i];
                aF[j][i] = scdiag * au[j][i] - scx * (uw + ue) - scy * (us + un) - darea * user->f_rhs(x, y, 0.0, user);
            }
        }
    }
    PetscCall(PetscLogFlops(11.0 * info->xm * info->ym));
    return 0;
}
// ENDFORM2DFUNCTION

PetscErrorCode Poisson3DFunctionLocal(DMDALocalInfo *info, PetscReal ***au,
                                      PetscReal ***aF, PoissonCtx *user)
{
    PetscInt i, j, k;
    PetscReal xyzmin[3], xyzmax[3], hx, hy, hz, dvol, scx, scy, scz, scdiag,
        x, y, z, ue, uw, un, us, uu, ud;
    PetscCall(DMGetBoundingBox(info->da, xyzmin, xyzmax));
    hx = (xyzmax[0] - xyzmin[0]) / (info->mx - 1);
    hy = (xyzmax[1] - xyzmin[1]) / (info->my - 1);
    hz = (xyzmax[2] - xyzmin[2]) / (info->mz - 1);
    dvol = hx * hy * hz;
    scx = user->cx * dvol / (hx * hx);
    scy = user->cy * dvol / (hy * hy);
    scz = user->cz * dvol / (hz * hz);
    scdiag = 2.0 * (scx + scy + scz);
    for (k = info->zs; k < info->zs + info->zm; k++)
    {
        z = xyzmin[2] + k * hz;
        for (j = info->ys; j < info->ys + info->ym; j++)
        {
            y = xyzmin[1] + j * hy;
            for (i = info->xs; i < info->xs + info->xm; i++)
            {
                x = xyzmin[0] + i * hx;
                if (i == 0 || i == info->mx - 1 || j == 0 || j == info->my - 1 || k == 0 || k == info->mz - 1)
                {
                    aF[k][j][i] = au[k][j][i] - user->g_bdry(x, y, z, user);
                    aF[k][j][i] *= scdiag;
                }
                else
                {
                    ue = (i + 1 == info->mx - 1) ? user->g_bdry(x + hx, y, z, user)
                                                 : au[k][j][i + 1];
                    uw = (i - 1 == 0) ? user->g_bdry(x - hx, y, z, user)
                                      : au[k][j][i - 1];
                    un = (j + 1 == info->my - 1) ? user->g_bdry(x, y + hy, z, user)
                                                 : au[k][j + 1][i];
                    us = (j - 1 == 0) ? user->g_bdry(x, y - hy, z, user)
                                      : au[k][j - 1][i];
                    uu = (k + 1 == info->mz - 1) ? user->g_bdry(x, y, z + hz, user)
                                                 : au[k + 1][j][i];
                    ud = (k - 1 == 0) ? user->g_bdry(x, y, z - hz, user)
                                      : au[k - 1][j][i];
                    aF[k][j][i] = scdiag * au[k][j][i] - scx * (uw + ue) - scy * (us + un) - scz * (uu + ud) - dvol * user->f_rhs(x, y, z, user);
                }
            }
        }
    }
    PetscCall(PetscLogFlops(14.0 * info->xm * info->ym * info->zm));
    return 0;
}

PetscErrorCode Poisson1DJacobianLocal(DMDALocalInfo *info, PetscScalar *au, Mat J, Mat Jpre, PoissonCtx *user)
{ // PetscScala用于表示实数或复数
    PetscInt i, ncols;
    PetscReal xmin[1], xmax[1], h, v[3];
    MatStencil col[3], row;

    PetscCall(DMGetBoundingBox(info->da, xmin, xmax));
    h = (xmax[0] - xmin[0]) / (info->mx - 1);
    for (i = info->xs; i < info->xs + info->xm; i++)
    {
        row.i = i;
        col[0].i = i;
        ncols = 1;
        if (i == 0 || i == info->mx - 1)
        {
            v[0] = user->cx * 2.0 / h;
        }
        else
        {
            v[0] = user->cx * 2.0 / h;
            if (i - 1 > 0)
            {
                col[ncols].i = i - 1;
                v[ncols++] = -user->cx / h;
            }
            if (i + 1 < info->mx - 1)
            {
                col[ncols].i = i + 1;
                v[ncols++] = -user->cx / h;
            }
        }
        PetscCall(MatSetValuesStencil(Jpre, 1, &row, ncols, col, v, INSERT_VALUES));
    }

    PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));
    if (J != Jpre)
    {
        PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));
    }
    return 0;
}

PetscErrorCode Poisson2DJacobianLocal(DMDALocalInfo *info, PetscScalar **au,
                                      Mat J, Mat Jpre, PoissonCtx *user)
{
    PetscReal xymin[2], xymax[2], hx, hy, scx, scy, scdiag, v[5];
    PetscInt i, j, ncols;
    MatStencil col[5], row;

    PetscCall(DMGetBoundingBox(info->da, xymin, xymax));
    hx = (xymax[0] - xymin[0]) / (info->mx - 1);
    hy = (xymax[1] - xymin[1]) / (info->my - 1);
    scx = user->cx * hy / hx;
    scy = user->cy * hx / hy;
    scdiag = 2.0 * (scx + scy);
    for (j = info->ys; j < info->ys + info->ym; j++)
    {
        row.j = j;
        col[0].j = j;
        for (i = info->xs; i < info->xs + info->xm; i++)
        {
            row.i = i;
            col[0].i = i;
            ncols = 1;
            v[0] = scdiag;
            if (i > 0 && i < info->mx - 1 && j > 0 && j < info->my - 1)
            {
                if (i - 1 > 0)
                {
                    col[ncols].j = j;
                    col[ncols].i = i - 1;
                    v[ncols++] = -scx;
                }
                if (i + 1 < info->mx - 1)
                {
                    col[ncols].j = j;
                    col[ncols].i = i + 1;
                    v[ncols++] = -scx;
                }
                if (j - 1 > 0)
                {
                    col[ncols].j = j - 1;
                    col[ncols].i = i;
                    v[ncols++] = -scy;
                }
                if (j + 1 < info->my - 1)
                {
                    col[ncols].j = j + 1;
                    col[ncols].i = i;
                    v[ncols++] = -scy;
                }
            }
            PetscCall(MatSetValuesStencil(Jpre, 1, &row, ncols, col, v, INSERT_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));
    if (J != Jpre)
    {
        PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));
    }
    return 0;
}

PetscErrorCode Poisson3DJacobianLocal(DMDALocalInfo *info, PetscScalar ***au,
                                      Mat J, Mat Jpre, PoissonCtx *user)
{
    PetscReal xyzmin[3], xyzmax[3], hx, hy, hz, dvol, scx, scy, scz, scdiag, v[7];
    PetscInt i, j, k, ncols;
    MatStencil col[7], row;

    PetscCall(DMGetBoundingBox(info->da, xyzmin, xyzmax));
    hx = (xyzmax[0] - xyzmin[0]) / (info->mx - 1);
    hy = (xyzmax[1] - xyzmin[1]) / (info->my - 1);
    hz = (xyzmax[2] - xyzmin[2]) / (info->mz - 1);
    dvol = hx * hy * hz;
    scx = user->cx * dvol / (hx * hx);
    scy = user->cy * dvol / (hy * hy);
    scz = user->cz * dvol / (hz * hz);
    scdiag = 2.0 * (scx + scy + scz);
    for (k = info->zs; k < info->zs + info->zm; k++)
    {
        row.k = k;
        col[0].k = k;
        for (j = info->ys; j < info->ys + info->ym; j++)
        {
            row.j = j;
            col[0].j = j;
            for (i = info->xs; i < info->xs + info->xm; i++)
            {
                row.i = i;
                col[0].i = i;
                ncols = 1;
                v[0] = scdiag;
                if (i > 0 && i < info->mx - 1 && j > 0 && j < info->my - 1 && k > 0 && k < info->mz - 1)
                {
                    if (i - 1 > 0)
                    {
                        col[ncols].k = k;
                        col[ncols].j = j;
                        col[ncols].i = i - 1;
                        v[ncols++] = -scx;
                    }
                    if (i + 1 < info->mx - 1)
                    {
                        col[ncols].k = k;
                        col[ncols].j = j;
                        col[ncols].i = i + 1;
                        v[ncols++] = -scx;
                    }
                    if (j - 1 > 0)
                    {
                        col[ncols].k = k;
                        col[ncols].j = j - 1;
                        col[ncols].i = i;
                        v[ncols++] = -scy;
                    }
                    if (j + 1 < info->my - 1)
                    {
                        col[ncols].k = k;
                        col[ncols].j = j + 1;
                        col[ncols].i = i;
                        v[ncols++] = -scy;
                    }
                    if (k - 1 > 0)
                    {
                        col[ncols].k = k - 1;
                        col[ncols].j = j;
                        col[ncols].i = i;
                        v[ncols++] = -scz;
                    }
                    if (k + 1 < info->mz - 1)
                    {
                        col[ncols].k = k + 1;
                        col[ncols].j = j;
                        col[ncols].i = i;
                        v[ncols++] = -scz;
                    }
                }
                PetscCall(MatSetValuesStencil(Jpre, 1, &row, ncols, col, v, INSERT_VALUES));
            }
        }
    }
    PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));
    if (J != Jpre)
    {
        PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));
    }
    return 0;
}

PetscErrorCode InitialState(DM da, InitialType it, PetscBool gbdry, Vec u, PoissonCtx *user)
{
    DMDALocalInfo info;
    PetscRandom rctx;
    switch (it)
    {
    case ZEROS:
        PetscCall(VecSet(u, 0.0));
        break;
    case RANDOM:
        PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rctx));
        PetscCall(VecSetRandom(u, rctx));
        PetscCall(PetscRandomDestroy(&rctx));
        break;
    default:
        SETERRQ(PETSC_COMM_SELF, 4, "invalid InitialType ... how did I get here?\n");
    }
    if (!gbdry)
    {
        return 0;
    }
    PetscCall(DMDAGetLocalInfo(da, &info));
    switch (info.dim)
    {
    case 1:
    {
        PetscInt i;
        PetscReal xmax[1], xmin[1], h, x, *au;
        PetscCall(DMDAVecGetArray(da, u, &au));
        PetscCall(DMGetBoundingBox(da, xmin, xmax));
        h = (xmax[0] - xmin[0]) / (info.mx - 1);
        for (i = info.xs; i < info.xs + info.xm; i++)
        {
            if (i == 0 || i == info.mx - 1)
            {
                x = xmin[0] + i * h;
                au[i] = user->g_bdry(x, 0.0, 0.0, user);
            }
        }
        PetscCall(DMDAVecRestoreArray(da, u, &au));
        break;
    }
    case 2:
    {
        PetscInt i, j;
        PetscReal xymin[2], xymax[2], hx, hy, x, y, **au;
        PetscCall(DMDAVecGetArray(da, u, &au));
        PetscCall(DMGetBoundingBox(da, xymin, xymax));
        hx = (xymax[0] - xymin[0]) / (info.mx - 1);
        hy = (xymax[1] - xymin[1]) / (info.my - 1);
        for (j = info.ys; j < info.ys + info.ym; j++)
        {
            y = xymin[1] + j * hy;
            for (i = info.xs; i < info.xs + info.xm; i++)
            {
                if (i == 0 || i == info.mx - 1 || j == 0 || j == info.my - 1)
                {
                    x = xymin[0] + i * hx;
                    au[j][i] = user->g_bdry(x, y, 0.0, user);
                }
            }
        }
        PetscCall(DMDAVecRestoreArray(da, u, &au));
        break;
    }
    case 3:
    {
        PetscInt i, j, k;
        PetscReal xyzmin[3], xyzmax[3], hx, hy, hz, x, y, z, ***au;
        PetscCall(DMDAVecGetArray(da, u, &au));
        PetscCall(DMGetBoundingBox(da, xyzmin, xyzmax));
        hx = (xyzmax[0] - xyzmin[0]) / (info.mx - 1);
        hy = (xyzmax[1] - xyzmin[1]) / (info.my - 1);
        hz = (xyzmax[2] - xyzmin[2]) / (info.mz - 1);
        for (k = info.zs; k < info.zs + info.zm; k++)
        {
            z = xyzmin[2] + k * hz;
            for (j = info.ys; j < info.ys + info.ym; j++)
            {
                y = xyzmin[1] + j * hy;
                for (i = info.xs; i < info.xs + info.xm; i++)
                {
                    if (i == 0 || i == info.mx - 1 || j == 0 || j == info.my - 1 || k == 0 || k == info.mz - 1)
                    {
                        x = xyzmin[0] + i * hx;
                        au[k][j][i] = user->g_bdry(x, y, z, user);
                    }
                }
            }
        }
        PetscCall(DMDAVecRestoreArray(da, u, &au));
        break;
    }
    default:
        SETERRQ(PETSC_COMM_SELF, 5, "invalid dim from DMDALocalInfo\n");
    }
    return 0;
}