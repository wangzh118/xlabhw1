cmake_minimum_required(VERSION 3.10)

# 设置工程名称
project(MatrixOperations)

# 设置 C 标准为 C99
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# 包含头文件目录
include_directories(inc/algebra.h)

# 添加源文件
set(SOURCES
    src/main.c
    src/algebra.c
)

# 生成可执行文件
add_executable(MatrixOperations ${SOURCES})