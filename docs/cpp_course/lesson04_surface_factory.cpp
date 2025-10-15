#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

// 继续沿用前面课程的 Surface 抽象，这里聚焦于“按配置字符串构建派生对象”。

class Surface {
public:
  virtual ~Surface() = default;
  virtual std::string description() const = 0;
};

class ZPlane final : public Surface {
public:
  explicit ZPlane(double z) : z_{z} {}
  std::string description() const override { return "ZPlane(z=" + std::to_string(z_) + ")"; }

private:
  double z_;
};

class Sphere final : public Surface {
public:
  Sphere(double r, double cx, double cy, double cz)
    : r_{r}, cx_{cx}, cy_{cy}, cz_{cz}
  {
  }

  std::string description() const override
  {
    return "Sphere(r=" + std::to_string(r_) + ", center=" + std::to_string(cx_) + "," +
           std::to_string(cy_) + "," + std::to_string(cz_) + ")";
  }

private:
  double r_;
  double cx_;
  double cy_;
  double cz_;
};

class SurfaceFactory {
public:
  using Creator = std::function<std::unique_ptr<Surface>(const std::map<std::string, double>&)>;

  void register_type(std::string type, Creator creator)
  {
    registry_.emplace(std::move(type), std::move(creator));
  }

  std::unique_ptr<Surface> create(const std::string& type,
                                  const std::map<std::string, double>& params) const
  {
    auto it = registry_.find(type);
    if (it == registry_.end()) {
      throw std::runtime_error("Unknown surface type: " + type);
    }
    return it->second(params);
  }

private:
  std::map<std::string, Creator> registry_;
};

int main()
{
  SurfaceFactory factory;
  factory.register_type("z-plane", [](const std::map<std::string, double>& params) {
    return std::make_unique<ZPlane>(params.at("z0"));
  });
  factory.register_type("sphere", [](const std::map<std::string, double>& params) {
    return std::make_unique<Sphere>(params.at("radius"), params.at("cx"), params.at("cy"),
                                    params.at("cz"));
  });

  const std::pair<std::string, std::map<std::string, double>> configs[] = {
    {"z-plane", {{"z0", 10.0}}},
    {"sphere", {{"radius", 3.0}, {"cx", 0.0}, {"cy", 0.0}, {"cz", 5.0}}},
  };

  for (const auto& [type, params] : configs) {
    auto surface = factory.create(type, params);
    std::cout << "Constructed " << surface->description() << '\n';
  }
}
