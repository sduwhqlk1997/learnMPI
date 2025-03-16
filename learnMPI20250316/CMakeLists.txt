cmake_minimum_required(VERSION 3.22)
project(mpi_test)
find_package(MPI REQUIRED)

add_executable(main main.cpp)
# 将MPI库链接到可执行文件
target_link_libraries(main PRIVATE MPI::MPI_CXX)