# 开始学习非结构网格

# 需要链接的文件
ln -s /root/Cppackage/petsc/lib/petsc/bin/petsc_conf.py

ln -s /root/Cppackage/petsc/lib/petsc/bin/PetscBinaryIO.py
# gmsh 按2.2版本的格式生成网格的命令
gmsh ./meshes/trap.geo -2 -format msh2 -o trap.msh
# 转化为PETSC的网格格式
python3 ./mesh2petsc.py output.msh