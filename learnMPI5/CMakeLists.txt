cmake_minimum_required(VERSION 3.22)
project(learnMPI1)
find_package(MPI REQUIRED) # REQUIRED 的作用是确保MPI库是构建过程中必须依赖的

add_executable(main main.cpp)

target_link_libraries(main PRIVATE MPI::MPI_CXX) 