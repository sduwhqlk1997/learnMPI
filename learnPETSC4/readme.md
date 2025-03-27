# 代码内容
线性方程组求解

# 笔记
`KSPSetOperators(ksp,A,A)` ：其中第一个 $A$ 定义了线性算子；第二个 $A$ 表示在运行时会使用 preconditioner 

都设置好后，使用 `KSPSolve(ksp,b,x)` 求解线性方程组