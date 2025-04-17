# 设置CMake的最低版本要求
cmake_minimum_required(VERSION 3.10)

# 设置项目名称和使用的C++标准
project(MyCppProject VERSION 1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 添加头文件路径
include_directories(head)

# 添加源文件路径
file(GLOB SOURCE_FILES "source/*.cpp")

# 添加主程序文件
set(MAIN_FILE main.cpp)

# 将所有源文件添加到可执行文件
add_executable(main ${MAIN_FILE} ${SOURCE_FILES})