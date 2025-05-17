# 代码内容
有限元 2D poisson 不做 matrix free
尝试怎加调整区域尺寸和网格大小的命令

# 成功运行命令
`mpiexec -n 2 ./main -da_refine 4 -ksp_type cg -pc_type bjacobi -sub_pc_type icc -ksp_converged_reason -ksp_rtol 1e-10`

# 经验教训
1. 同一个Mat或者Vec只能进行修改或者读取，两个操作不能同步进行，若要同步进行则需要创建一个副本
2. 使用 `DMDAVecGetArray` 创建的局部向量副本不包含虚拟点，若要对虚拟点处的值修改则需要进程间通讯，如使用 `VecSetValues(b,NbInsert,indexbInsert,valuebInsert,INSERT_VALUES)`