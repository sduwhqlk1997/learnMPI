static char help[] = "A structured-grid Possion solver using DMDA+KSP.\n.\n";

# include <petsc.h>

extern PetscErrorCode formMatrix(DM, Mat); // 组装矩阵，PetscErrorCode数据类型，是一个整数类型，返回0表示成功，非0值表示错误
extern PetscErrorCode formExact(DM, Vec); 
extern PetscErrorCode formRHS(DM,Vec);

int main(int argc, char* argv[]){

}

// function definition
PetscErrorCode formMatrix(DM da, Mat A){
    DMDALocalInfo info; // 记录网格信息
    MatStencil row, col[5];
    PetscReal hx,hy,v[5];
    PetscInt i,j,ncols;

    PetscCall(DMDAGetLocalInfo(da,&info)); // 获取 `da` 中的网格信息
    hx = 1.0/(info.mx-1); hy = 1.0/(info.my-1); // xm,ym,zx: the size of the local mesh of the current process
    for (j = info.ys;j<info.ys+info.ym;j++){ // xs,ys,zs: the start index of the current process
        for(i=info.xs;i<info.xs+info.xm;i++){
            row.j = j; // 矩阵 A 关于点 (x_i,y_j) 的行
            row.i = i;
            col[0].j = j; // diagonal entry
            col[0].i = i;
            ncols = 1;
            if (i==0 || i==info.mx-1 || j==0 || j==info.my-1){
                v[0] = 1.0; // 边界上的方程
            }else {
                v[0] = 2*(hy/hx + hx/hy); // 对角线上的值？
                if (i-1>0){
                    col[ncols].j=j; col[ncols].i = i-1;
                    v[ncols++]=-hy/hx;
                }
                if (i+1<info.mx-1){
                    col[ncols].j = j; col[ncols].i = i+1;
                    v[ncols++] = -hy/hx;
                }
                if(j-1>0){
                    col[ncols].j = j-1; col[ncols].i = i;
                    v[ncols++] = -hx/hy;
                }
                if(j+1 < info.my-1) {
                    col[ncols].j = j+1; col[ncols].i = i;
                    v[ncols++] = -hx/hy;
                }
            }
            PetscCall(MatSetValuesStencil(A,1,&row,ncols,col,v,INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
    return 0;
}

PetscErrorCode formExact(DM da, Vec uexact){
    PetscInt i,j;
    PetscReal hx, hy, x, y, **auexact;
    DMDALocalInfo info;

    PetscCall(DMDAGetLocalInfo(da,&info));
    hx = 1.0/(info.mx-1); hy = 1.0/(info.my-1);
    PetscCall(DMDAVecGetArray(da, uexact, &auexact));
    for (j = info.ys; j< info.ys+info.ym; j++){
        y = j * hy;
        for(i = info.xs; i<info.xs+info.xm;i++){
            x = i*hx;
            auexact[j][i]=x*x*(1.0-x*x) * y*y * (y*y-1.0);
        }
    }
    PetscCall(DMDAVecRestoreArray(da,uexact,&auexact));
    return 0;
}