#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

// 本课强调“多态材料 + 缓存”理念：
// - MaterialModel 抽象出求宏观截面的接口。
// - ConstantMaterial 与 TemperatureDependentMaterial 展示策略模式。
// - CachedMaterialWrapper 重用上次能量结果，仿照 OpenMC 的 xs 缓存。

class MaterialModel {
public:
  virtual ~MaterialModel() = default;
  virtual double macroscopic_xs(double energy_eV) const = 0;
  virtual std::string name() const = 0;
};

class ConstantMaterial final : public MaterialModel {
public:
  ConstantMaterial(std::string id, double xs) : id_{std::move(id)}, xs_{xs} {}

  double macroscopic_xs(double) const override { return xs_; }
  std::string name() const override { return id_; }

private:
  std::string id_;
  double xs_;
};

class TemperatureDependentMaterial final : public MaterialModel {
public:
  TemperatureDependentMaterial(std::string id, double base_xs, double temperature)
    : id_{std::move(id)}
    , base_xs_{base_xs}
    , temperature_{temperature}
  {
  }

  double macroscopic_xs(double energy_eV) const override
  {
    double shift = std::sqrt(energy_eV / 1e5);
    return base_xs_ * (1.0 + 0.1 * std::sin(shift)) * (1.0 + 0.001 * temperature_);
  }

  std::string name() const override { return id_; }

private:
  std::string id_;
  double base_xs_;
  double temperature_;
};

class CachedMaterialWrapper final : public MaterialModel {
public:
  explicit CachedMaterialWrapper(std::shared_ptr<MaterialModel> impl)
    : impl_{std::move(impl)}
  {
  }

  double macroscopic_xs(double energy_eV) const override
  {
    if (!last_energy_.has_value() || std::abs(*last_energy_ - energy_eV) > 1e-6) {
      last_energy_ = energy_eV;
      last_value_ = impl_->macroscopic_xs(energy_eV);
      cache_hits_ = 0;
    } else {
      ++cache_hits_;
    }
    return *last_value_;
  }

  std::string name() const override { return impl_->name(); }

  int cache_hits() const { return cache_hits_; }

private:
  std::shared_ptr<MaterialModel> impl_;
  mutable std::optional<double> last_energy_;
  mutable std::optional<double> last_value_;
  mutable int cache_hits_ = 0;
};

int main()
{
  auto fuel = std::make_shared<CachedMaterialWrapper>(
    std::make_shared<TemperatureDependentMaterial>("UO2", 0.6, 900.0));
  auto water = std::make_shared<CachedMaterialWrapper>(
    std::make_shared<ConstantMaterial>("H2O", 0.2));

  std::vector<double> energies{1e3, 2e3, 2e3, 5e3};
  for (double E : energies) {
    std::cout << "fuel xs at " << E << " eV = " << fuel->macroscopic_xs(E)
              << " (cache hits: " << fuel->cache_hits() << ")\n";
    std::cout << "water xs at " << E << " eV = " << water->macroscopic_xs(E)
              << " (cache hits: " << water->cache_hits() << ")\n";
  }
}
