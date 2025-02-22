cmake_minimum_required(VERSION 3.10)
project(fmt_test)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找 fmt 库
find_package(fmt REQUIRED)

# 指定源文件
add_executable(fmt_test test_fmt.cpp)

# 链接 fmt 库
target_link_libraries(fmt_test fmt::fmt)
