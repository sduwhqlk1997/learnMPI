# 代码内容
线性方程组求解

# 笔记
1. `KSPSetOperators(ksp,A,A)` ：其中第一个 $A$ 定义了线性算子；第二个 $A$ 表示在运行时会使用 preconditioner 
2. 都设置好后，使用 `KSPSolve(ksp,b,x)` 求解线性方程组
3. `PetscCall(VecSetValues(b,4,j,ab,INSERT_VALUES));` 中的`INSERT_VALUES`表示将提供的值直接插入到向量的指定位置，覆盖该位置原有的值。