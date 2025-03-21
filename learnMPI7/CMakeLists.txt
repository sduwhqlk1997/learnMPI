cmake_minimum_required(VERSION 3.22)
project(learnMPI1)
set(CMAKE_BUILD_TYPE Debug) # 生成debug符号
add_definitions(-DDEBUG) # 定义 DEBUG 宏

find_package(MPI REQUIRED) # REQUIRED 的作用是确保MPI库是构建过程中必须依赖的

add_executable(main main.cpp)

target_link_libraries(main PRIVATE MPI::MPI_CXX) 