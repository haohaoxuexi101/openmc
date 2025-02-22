#include "lattice.h"
#include "pugixml.hpp" // 需要引入XML解析库
#include <cassert>
#include <iostream>

int main()
{
  try {
    // 加载XML节点
    pugi::xml_document doc;
    doc.load_string(R"(
            <lattice>
                <id>1</id>
                <dimension>3 3</dimension>
                <lower_left>0.0 0.0</lower_left>
                <pitch>1.0 1.0</pitch>
                <universes>1 1 1 1 1 1 1 1 1</universes>
            </lattice>
        )");

    pugi::xml_node lat_node = doc.child("lattice");

    // 创建一个2D矩形晶格
    std::cout << "Creating a 2D rectangular lattice..." << std::endl;
    openmc::RectLattice lattice(lat_node);

    // 设置晶格的参数并调用distance函数进行测试
    std::cout << "Lattice created." << std::endl;

    // 定义粒子的位置和方向
    openmc::Position r = {0.5, 0.5, 0.0}; // 粒子位于第一个晶格中心
    // openmc::Position r = {1.0, 0.5, 0.0}; // 直接位于晶格的右边界附近
    openmc::Direction u = {1.0, 1.0, 0.0}; // 粒子沿x方向移动

    // 当前晶格坐标
    openmc::array<int, 3> i_xyz = {1, 1, 0}; // 位于晶格(1,1)位置

    // 调用distance函数计算距离和方向
    std::cout << "Calculating distance to next boundary..." << std::endl;
    auto result = lattice.distance(r, u, i_xyz);

    // 打印结果
    std::cout << "Distance to next boundary: " << result.first << std::endl;
    std::cout << "Lattice transition direction: {" << result.second[0] << ", "
              << result.second[1] << ", " << result.second[2] << "}"
              << std::endl;

    // 验证结果是否正确
    assert(result.first > 0.0 && "Distance should be positive");
    // assert(result.second[0] == 1 && result.second[1] == 0 && "Lattice should
    // move in x direction");

    std::cout << "Test passed!" << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "An error occurred during testing: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
