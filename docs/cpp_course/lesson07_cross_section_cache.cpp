#include <cmath>
#include <iostream>
#include <map>
#include <tuple>

struct MaterialXsKey {
  int material_id;
  int temperature_id;
  double energy_group;

  bool operator<(const MaterialXsKey& other) const
  {
    return std::tie(material_id, temperature_id, energy_group) <
           std::tie(other.material_id, other.temperature_id, other.energy_group);
  }
};

class CrossSectionLibrary {
public:
  double query(int material, int temperature, double energy_group) const
  {
    return 0.2 * material + 0.05 * temperature + std::log1p(energy_group);
  }
};

class CrossSectionCache {
public:
  explicit CrossSectionCache(std::size_t capacity) : capacity_{capacity} {}

  double get_or_compute(const MaterialXsKey& key, const CrossSectionLibrary& lib)
  {
    if (auto it = cache_.find(key); it != cache_.end()) {
      std::cout << "cache hit\n";
      return it->second;
    }
    if (cache_.size() == capacity_) {
      cache_.erase(cache_.begin());
      std::cout << "evict oldest entry\n";
    }
    const double xs = lib.query(key.material_id, key.temperature_id, key.energy_group);
    cache_.emplace(key, xs);
    std::cout << "cache miss -> insert\n";
    return xs;
  }

private:
  std::size_t capacity_;
  std::map<MaterialXsKey, double> cache_;
};

int main()
{
  CrossSectionLibrary library;
  CrossSectionCache cache{2};

  MaterialXsKey keys[] = {
    {1, 600, 0.5},
    {1, 600, 0.5},
    {2, 900, 1.0},
    {3, 1200, 1.5},
    {1, 600, 0.5},
  };

  for (const auto& key : keys) {
    double xs = cache.get_or_compute(key, library);
    std::cout << "  xs=" << xs << '\n';
  }
}
