# 代码内容
三对角线性方程组
## 程序运行命令
`$time mpiexec -n P ./main -tri_m 20000000 -ksp_rtol 1.0e-10 -ksp_type X -pc_type Y` 这里 `X` 和 `Y` 是待定的迭代法类型和预条件，可选择如下组合：
1. -ksp_type preonly -pc_type lu
2. -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu
3. -ksp_type preonly -pc_type cholesky
4. -ksp_type cg -pc_type icc
   
如：time mpiexec -n 8 ./main -tri_m 20000000 -ksp_rtol 1.0e-10 -ksp_type cg -pc_type bjacobi -sub_pc_type icc
# 笔记
1. `MatSetSizes (A,PETSC_DECIDE,PETSC_DECIDE,m,m)` ：第2、3个变量的`PETSC_DECIDE`为局部大小，会根据全局大小和进程数来决定；最后两个变量`m`为全局大小，表示矩阵$A$的为$m\times m$的方阵。
2. `MatGetOwnershipRange (A,&Istart ,&Iend )`：用于获取当前进程的矩阵A的首行`Istart`和末行`Iend`。可用来并行组装矩阵
3. 在`PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY))`中 `MAT_FINAL_ASSEMBLY` 标志矩阵的最终组装状态；矩阵 $A$ 将被标记为“已组装”，并且其内部数据结构将被优化，以便后续操作可以高效执行。在矩阵组装过程中，PETSc可能会使用一些临时数据结构来存储未组装的元素。`MAT_FINAL_ASSEMBLY` 会告诉 PETSc 释放这些临时数据结构，以节省内存。
