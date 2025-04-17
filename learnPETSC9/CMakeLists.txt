cmake_minimum_required(VERSION 3.10)
project(PETScTest)
# Debug
#set(CMAKE_BUILD_TYPE Debug) # 生成debug符号
#add_definitions(-DDEBUG) # 定义 DEBUG 宏
# 设置 PETSc 的安装路径
set(PETSc_DIR "/root/Cppackage/petsc")  # 替换为你的 PETSc 安装路径

#头文件路径
include_directories(${PETSc_DIR}/include head)
# 添加源文件路径
file(GLOB SOURCE_FILES "source/*.cpp")
# 添加可执行文件
add_executable(main main.cpp ${SOURCE_FILES})

# 查找MPI库
find_package(MPI REQUIRED)

# 包含 MPI 头文件路径
target_include_directories(main PRIVATE ${MPI_INCLUDE_PATH})

# 链接 PETSc 库
target_link_directories(main PRIVATE ${PETSc_DIR}/lib)
target_link_libraries(main PRIVATE petsc ${MPI_C_LIBRARIES})