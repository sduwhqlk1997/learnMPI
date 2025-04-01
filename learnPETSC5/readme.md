# 代码内容
三对角线性方程组
# 笔记
1. `MatSetSizes (A,PETSC_DECIDE,PETSC_DECIDE,m,m)` ：第2、3个变量的`PETSC_DECIDE`为局部大小，会根据全局大小和进程数来决定；最后两个变量`m`为全局大小，表示矩阵$A$的为$m\times m$的方阵。
2. `MatGetOwnershipRange (A,&Istart ,&Iend )`：用于获取当前进程的矩阵A的首行`Istart`和末行`Iend`。