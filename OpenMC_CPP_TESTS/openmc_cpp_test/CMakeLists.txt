cmake_minimum_required(VERSION 3.10)

# 项目名称
project(OpenMC_Test)

# 指定 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 查找 HDF5 库
find_package(HDF5 REQUIRED COMPONENTS CXX)

# 定义多个源文件名
set(SOURCE_FILES
  test_particle
  test_lattice
  # 在这里添加其他源文件名
)

# 遍历源文件列表，生成多个可执行文件
foreach(file ${SOURCE_FILES})
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} PRIVATE openmc ${HDF5_LIBRARIES})
endforeach()

# 指定包含目录，以便找到 HDF5 和 OpenMC 的头文件
include_directories(${HDF5_INCLUDE_DIRS} /home/ssn/ssn_mc/openmc/include/openmc)  # 替换为 OpenMC 头文件的实际路径
# include_directories(${HDF5_INCLUDE_DIRS})  # 替换为 OpenMC 头文件的实际路径
