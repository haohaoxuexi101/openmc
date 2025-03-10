cmake_minimum_required(VERSION 3.10)

# 项目名称
project(OpenMC_Test)

# 指定 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 查找 HDF5 库
find_package(HDF5 REQUIRED COMPONENTS CXX)

# 查找所有的 .cpp 文件
file(GLOB SOURCES "*.cpp")

# 如果你想要在子目录中也进行递归查找 .cpp 文件，可以使用下面这行
# file(GLOB_RECURSE SOURCES "src/*.cpp")

# 使用 link_directories 指定库路径
# 与target_link_libraries中的openmc对应，为其提供路径
link_directories(/home/ssn/ssn_mc/openmc/build/lib)

# 遍历每个源文件，生成目标文件
foreach(file ${SOURCES})
  # 获取文件名（不带扩展名）
  get_filename_component(target_name ${file} NAME_WE)

  # 生成可执行文件
  add_executable(${target_name} ${file})

  # # 链接 OpenMC 和 HDF5 库
  # # 与link_directories中的路径对应，在其路径中查找库
  # # libopenmc.so或者openmc都可以
  # # CMake 会自动为库名添加 lib 前缀和 .so 后缀
  # target_link_libraries(${target_name} PRIVATE libopenmc.so ${HDF5_LIBRARIES})

  # 链接 OpenMC 和 HDF5 库
  target_link_libraries(${target_name} PRIVATE /home/ssn/ssn_mc/openmc/build/lib/libopenmc.so ${HDF5_LIBRARIES})

  # 为目标设置包含目录
  target_include_directories(${target_name} PRIVATE
    # 由于所有头文件都是openmc/
    # 因此该路径只取openmc上一层，只到openmc
    /home/ssn/ssn_mc/openmc/include  # OpenMC 头文件路径
    # 非常重要，原来此CMakeLists.txt能够成功
    # 是因为/usr/local/include中含有对应的头文件
    /home/ssn/ssn_mc/openmc/build/include # 需要 version.h
    ${HDF5_INCLUDE_DIRS}  # HDF5 头文件路径
  )
endforeach()

# # 指定包含目录，以便找到 HDF5 和 OpenMC 的头文件
# include_directories(${HDF5_INCLUDE_DIRS} /home/ssn/ssn_mc/openmc/include/openmc)