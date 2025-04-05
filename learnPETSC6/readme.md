# 代码内容
求解二维Possion方程
## 程序运行命令
### 正常运行
`./main -da_refine 4 -ksp_monitor -ksp_monitor_solution draw -draw_pause 0.1`
### Debug

# 笔记
1. `DMDAVecGetArray(da, uexact, &auexact)` 和 `DMDAVecRestoreArray(da,uexact,&auexact)` 成对使用
2. DMDA 是 Distributed Array (DMDA) 的缩写，DMDA 是 PETSc 的 DM（Data Management）模块的一部分，专门用于处理与网格相关的数据，特别是在求解偏微分方程（PDE）时。
3. `DMDALocalInfo info`数据类型：
   1. `info.mx`, `info.my`：全局网格数 (GLOBAL)
   2. `info.xm`, `info.ym`：当前线程的局部网格数 (LOCAL)
   3. `info.xs`, `info.ys`：当前线程的网格起点 (LOCAL)
