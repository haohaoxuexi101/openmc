#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

// 课程 3：导航上下文（NavigationContext）的接口分离
// 借鉴说明：该课映射 OpenMC `src/geometry/navigator.h` 与
// `src/geometry/geometry.h` 中 Navigator/GeometryState 的协作模式，
// 特别强调定位调用与上下文记录的解耦；示例使用简化的字符串轨迹
// 来说明接口职责，而非直接复刻源码实现细节。

class Navigator;  // 前向声明，体现解耦

class NavigationContext {
public:
  NavigationContext(const std::string& name, std::shared_ptr<Navigator> nav)
      : name_(name), navigator_(std::move(nav)) {}

  void push_position(double x, double y, double z);

  void show_history() const {
    std::cout << "上下文 " << name_ << " 记录的定位轨迹:\n";
    for (const auto& item : history_) {
      std::cout << "  - " << item << "\n";
    }
  }

private:
  std::string name_;
  std::shared_ptr<Navigator> navigator_;
  std::vector<std::string> history_;
};

class Navigator {
public:
  std::string locate(double x, double y, double z) const {
    if (x * x + y * y + z * z < 0.25) return "fuel";
    if (std::abs(x) < 1.0 && std::abs(y) < 1.0 && std::abs(z) < 1.0)
      return "moderator";
    return "outside";
  }
};

void NavigationContext::push_position(double x, double y, double z) {
  std::string region = navigator_->locate(x, y, z);
  history_.push_back("位置(" + std::to_string(x) + ", " + std::to_string(y) +
                     ", " + std::to_string(z) + ") -> " + region);
}

int main() {
  auto navigator = std::make_shared<Navigator>();
  NavigationContext context("主几何", navigator);

  context.push_position(0.0, 0.0, 0.0);
  context.push_position(0.7, 0.1, 0.1);
  context.push_position(1.2, 0.0, 0.0);

  context.show_history();
  return 0;
}
