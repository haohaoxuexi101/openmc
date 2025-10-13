#include <cmath>
#include <iostream>

// 本示例聚焦 OpenMC 代码中最基础的 Vec3 工具类，演示如何封装三维向量算子。
// 在 OpenMC 主体中，粒子的位置、方向均依赖 Vec3。我们通过中文注释细致解释
// 每个成员函数的设计考虑，帮助读者掌握面向对象封装细节。

namespace lesson_vec3 {

// ============================= Vec3 向量类 =============================
// 该类与 OpenMC C++ 核心中的 Vec3 结构体同源：
//   1. 公开成员以便直接访问，保持性能；
//   2. 同时提供运算符，保证外部调用表达式自然。
struct Vec3 {
  double x = 0.0;  // x 分量
  double y = 0.0;  // y 分量
  double z = 0.0;  // z 分量

  Vec3() = default;  // 默认构造器保持零向量
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  // 运算符重载展示“值语义”——让向量加减乘法像数学式一样可读。
  Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
  Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

  // += 返回引用用于链式调用，和标准库容器保持一致接口风格。
  Vec3& operator+=(const Vec3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
};

// 工具函数使用独立的自由函数形式，体现“算法与数据分离”的设计理念。
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalise(const Vec3& v) {
  double n = norm(v);
  if (n == 0.0) {
    // 以 OpenMC 的防御式编程为例：当输入为零向量时，给出默认方向避免除零。
    return {0.0, 0.0, 1.0};
  }
  return {v.x / n, v.y / n, v.z / n};
}

}  // namespace lesson_vec3

int main() {
  using lesson_vec3::Vec3;
  using lesson_vec3::dot;
  using lesson_vec3::norm;
  using lesson_vec3::normalise;

  // ============================= 示例用法 =============================
  // 1. 构造向量并进行加法，体现值语义的便利。
  Vec3 position{1.0, 2.0, 3.0};
  Vec3 shift{0.0, -1.0, 0.5};
  Vec3 new_position = position + shift;

  // 2. 计算方向余弦，展示 normalise 与 dot 的组合使用。
  Vec3 direction = normalise({2.0, -1.0, 0.5});
  double cosine = dot(direction, {0.0, 0.0, 1.0});

  // 3. 输出结果，验证接口直观性。
  std::cout << "归一化方向 = (" << direction.x << ", " << direction.y << ", " << direction.z
            << ")\n";
  std::cout << "与 z 轴余弦 = " << cosine << "\n";
  std::cout << "平移后位置 = (" << new_position.x << ", " << new_position.y << ", "
            << new_position.z << ")\n";
  std::cout << "平移长度 = " << norm(shift) << "\n";
}
