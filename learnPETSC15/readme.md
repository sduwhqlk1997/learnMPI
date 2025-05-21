![alt text](image.png)
运行命令：
mpiexec -n 2 ./main -ksp_type cg -da_refine 6 -ksp_converged_reason -ksp_rtol 1e-10 -pc_type asm -pc_asm_overlap 4 -sub_pc_type icc
# 代码内容
有限元 

# 局部向量->全局向量
![alt text](image.png)
![alt text](image-1.png)
![alt text](image-2.png)