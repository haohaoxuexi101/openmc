#define private public  // 将 private 替换为 public
#include "particle.h"
#include "constants.h"
#include <iostream>
#include <cassert>

int main() {
    openmc::Particle particle;

    // 初始化粒子类型和能量
    particle.type() = openmc::ParticleType::neutron;
    particle.E() = 1.0e6;  // 1 MeV

    // 检查粒子速度是否正常
    double particle_speed = particle.speed();
    std::cout << "Particle speed (neutron, 1 MeV): " << particle_speed << " m/s" << std::endl;
    assert(particle_speed > 0.0 && "Speed should be positive for neutron");

    // 初始化坐标系统：确保coord_数组有足够的空间
    particle.coord_.resize(1);  // 确保coord_数组有一个元素

    particle.n_coord() = 1;  // 设置坐标层级
    particle.coord(0).reset();  // 重置并确保第一个坐标层级初始化

    // 检查坐标是否被初始化
    if (particle.n_coord() > 0) {
        std::cout << "Particle coordinate levels initialized: " << particle.n_coord() << std::endl;

        // 初始化坐标
        // particle.coord(0).r.x = 0.0;
        // particle.coord(0).r.y = 0.0;
        // particle.coord(0).r.z = 0.0;  // 设置初始位置
        particle.coord(0).r = {0, 0, 0};

        particle.coord(0).u.x = 1.0;
        particle.coord(0).u.y = 0.0;
        particle.coord(0).u.z = 0.0;  // 设置方向（沿 x 轴）
    } else {
        std::cerr << "Error: Particle coordinate not properly initialized!" << std::endl;
        return -1;
    }

    // 打印初始位置
    std::cout << "Initial particle position: ("
              << particle.coord(0).r.x << ", "
              << particle.coord(0).r.y << ", "
              << particle.coord(0).r.z << ")" << std::endl;

    // 移动粒子一段距离
    particle.move_distance(1.0);  // 移动 1 米
    std::cout << "Particle position after moving 1 meter: ("
              << particle.coord(0).r.x << ", "
              << particle.coord(0).r.y << ", "
              << particle.coord(0).r.z << ")" << std::endl;

    // 检查粒子是否在 x 轴上移动
    assert(particle.coord(0).r.x == 1.0 && "Particle should have moved 1 meter along x-axis");

    return 0;
}
