cmake_minimum_required(VERSION 3.22)
project(learnMPI1)
find_package(MPI REQUIRED)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE MPI::MPI_CXX)