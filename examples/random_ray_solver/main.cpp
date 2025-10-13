#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// =============================================================================
//  三维随机射线多群多区域求解器 (Random Ray Monte Carlo Transport Solver)
// -----------------------------------------------------------------------------
//  本程序源于 OpenMC 的面向对象设计理念，将材料(Material)、几何(Geometry)、
//  粒子(Particle)、随机源(Source)、计分(Tally) 与求解驱动(Solver) 解耦，实现
//  面向对象的组合式建模。程序支持：
//    - 多群截面 (任意群数)；
//    - 多区域三维轴对齐盒体几何 (嵌套/多层)；
//    - 体积分布的随机射线抽样与曲线长度计分；
//    - 吸收、散射、裂变的碰撞抽样与二次粒子生成；
//    - 可配置的粒子历史数、随机种子、最大追踪步数；
//    - 通过外部输入文件配置材料、几何与源项。 
//  该求解器定位为可以直接运行的教学/验证代码，帮助读者在掌握 OpenMC 设计风格
//  的同时，快速搭建多群多区域的蒙特卡洛输运实验。全文件提供中文注释，详细解释
//  每个类的职责与交互。若希望先逐模块学习，请参照同目录 `lessons/` 子文件夹中
//  的 6 个课程程序，那里分别讲解域构建、射线抽样、区域穿越、碰撞物理、裂变银行
//  与求解器驱动，随后再阅读本文件。
// =============================================================================

namespace random_ray_solver {

// ============================= 数学与向量工具 ================================

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  Vec3() = default;
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
  Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

  Vec3& operator+=(const Vec3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalize(const Vec3& v) {
  double n = norm(v);
  if (n == 0.0) return {0.0, 0.0, 1.0};
  return {v.x / n, v.y / n, v.z / n};
}

// ============================= 随机引擎封装 ==================================

class RandomEngine {
 public:
  explicit RandomEngine(unsigned seed = 20240602u) : rng_(seed) {}

  double uniform() { return uniform_(rng_); }

  Vec3 isotropic_direction() {
    double mu = 2.0 * uniform() - 1.0;
    double phi = 2.0 * kPi * uniform();
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), mu};
  }

  std::size_t pick(const std::vector<double>& pdf) {
    std::discrete_distribution<std::size_t> dist(pdf.begin(), pdf.end());
    return dist(rng_);
  }

 private:
  static constexpr double kPi = 3.14159265358979323846;
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

// ============================= 材料系统 ======================================

class Material {
 public:
  Material() = default;

  Material(std::string name, std::vector<double> sigma_a, std::vector<double> sigma_f,
           std::vector<double> nu, std::vector<std::vector<double>> sigma_s,
           std::vector<double> chi)
      : name_(std::move(name)), sigma_a_(std::move(sigma_a)), sigma_f_(std::move(sigma_f)),
        nu_(std::move(nu)), sigma_s_(std::move(sigma_s)), chi_(std::move(chi)) {
    std::size_t g = sigma_a_.size();
    if (sigma_f_.size() != g || nu_.size() != g || sigma_s_.size() != g || chi_.size() != g) {
      throw std::runtime_error("材料截面数组长度不一致");
    }
    for (const auto& row : sigma_s_) {
      if (row.size() != g) throw std::runtime_error("散射矩阵必须是方阵");
    }
  }

  [[nodiscard]] std::size_t groups() const noexcept { return sigma_a_.size(); }
  [[nodiscard]] const std::string& name() const noexcept { return name_; }

  [[nodiscard]] double sigma_a(std::size_t g) const { return sigma_a_.at(g); }
  [[nodiscard]] double sigma_f(std::size_t g) const { return sigma_f_.at(g); }
  [[nodiscard]] double nu(std::size_t g) const { return nu_.at(g); }
  [[nodiscard]] const std::vector<double>& chi() const noexcept { return chi_; }
  [[nodiscard]] const std::vector<double>& scatter_row(std::size_t g) const {
    return sigma_s_.at(g);
  }

  [[nodiscard]] double sigma_t(std::size_t g) const {
    double scatter = 0.0;
    for (double val : sigma_s_[g]) scatter += val;
    return sigma_a_[g] + sigma_f_[g] + scatter;
  }

 private:
  std::string name_;
  std::vector<double> sigma_a_;
  std::vector<double> sigma_f_;
  std::vector<double> nu_;
  std::vector<std::vector<double>> sigma_s_;
  std::vector<double> chi_;
};

// ============================= 几何系统 ======================================

struct AxisAlignedBox {
  Vec3 min;
  Vec3 max;

  [[nodiscard]] bool contains(const Vec3& p) const {
    return (p.x >= min.x && p.x <= max.x) && (p.y >= min.y && p.y <= max.y) &&
           (p.z >= min.z && p.z <= max.z);
  }

  [[nodiscard]] double distance_to_boundary(const Vec3& pos, const Vec3& dir) const {
    const double kInf = 1e300;
    double t_min = kInf;

    auto distance_axis = [&](double p, double d, double min_val, double max_val) {
      if (std::fabs(d) < 1e-12) return kInf;
      double t1 = (min_val - p) / d;
      double t2 = (max_val - p) / d;
      double t_near = std::min(t1, t2);
      double t_far = std::max(t1, t2);
      if (t_far <= 0.0) return kInf;
      if (t_near > 1e-12) return t_near;
      if (t_far > 1e-12) return t_far;
      return kInf;
    };

    t_min = std::min({distance_axis(pos.x, dir.x, min.x, max.x),
                      distance_axis(pos.y, dir.y, min.y, max.y),
                      distance_axis(pos.z, dir.z, min.z, max.z)});
    return t_min;
  }
};

struct Region {
  std::string name;
  AxisAlignedBox box;
  const Material* material = nullptr;
};

class Geometry {
 public:
  void add_region(Region region) { regions_.push_back(std::move(region)); }

  [[nodiscard]] const Region* locate(const Vec3& p) const {
    for (auto it = regions_.rbegin(); it != regions_.rend(); ++it) {
      if (it->box.contains(p)) return &(*it);
    }
    return nullptr;
  }

  [[nodiscard]] double distance_to_boundary(const Vec3& pos, const Vec3& dir,
                                            const Region& region) const {
    return region.box.distance_to_boundary(pos, dir);
  }

  [[nodiscard]] std::size_t region_index(const Region* region) const {
    for (std::size_t i = 0; i < regions_.size(); ++i) {
      if (&regions_[i] == region) return i;
    }
    throw std::runtime_error("未找到区域索引");
  }

  [[nodiscard]] const std::vector<Region>& regions() const noexcept { return regions_; }

 private:
  std::vector<Region> regions_;
};

// ============================= 粒子与粒子银行 =================================

struct Particle {
  Vec3 position;
  Vec3 direction;
  std::size_t group = 0;
  bool alive = true;
};

struct ParticleBank {
  std::vector<Particle> stack;

  void push(const Particle& p) { stack.push_back(p); }

  bool empty() const noexcept { return stack.empty(); }

  Particle pop() {
    Particle p = stack.back();
    stack.pop_back();
    return p;
  }
};

// ============================= Tallies 与统计 =================================

struct Tallies {
  std::vector<std::vector<double>> flux;         // 区域-能群曲线长度通量
  std::vector<std::vector<double>> absorption;   // 区域-能群吸收率
  std::vector<std::vector<double>> scattering;   // 区域-能群散射率
  std::vector<std::vector<double>> fission;      // 区域-能群裂变率
  double leakage = 0.0;                          // 漏出总权重
  double starting_histories = 0.0;               // 初始粒子历史数
};

// ============================= 配置数据结构 ==================================

struct MaterialInput {
  std::string name;
  std::vector<double> sigma_a;
  std::vector<double> sigma_f;
  std::vector<double> nu;
  std::vector<std::vector<double>> sigma_s;
  std::vector<double> chi;
};

struct RegionInput {
  std::string name;
  Vec3 min;
  Vec3 max;
  std::string material_name;
};

struct SourceInput {
  std::string region_name;              // 粒子抽样所在区域
  std::vector<double> group_pdf;        // 能群概率密度
};

struct GlobalOptions {
  std::size_t histories = 10000;        // 初始粒子数
  unsigned seed = 20240602u;
  std::size_t max_steps = 2000;         // 每个粒子的最大追踪步数
};

struct Config {
  std::size_t groups = 0;
  GlobalOptions options;
  std::vector<MaterialInput> materials;
  std::vector<RegionInput> regions;
  SourceInput source;
};

// ============================= 工具函数：字符串处理 ===========================

inline std::string trim(const std::string& s) {
  std::size_t first = 0;
  while (first < s.size() && std::isspace(static_cast<unsigned char>(s[first]))) ++first;
  std::size_t last = s.size();
  while (last > first && std::isspace(static_cast<unsigned char>(s[last - 1]))) --last;
  return s.substr(first, last - first);
}

inline std::vector<std::string> split(const std::string& s) {
  std::istringstream iss(s);
  std::vector<std::string> tokens;
  std::string token;
  while (iss >> token) tokens.push_back(token);
  return tokens;
}

inline std::vector<double> parse_double_list(const std::string& value) {
  std::vector<double> result;
  std::istringstream iss(value);
  double v = 0.0;
  while (iss >> v) result.push_back(v);
  return result;
}

inline Vec3 parse_vec3(const std::vector<std::string>& tokens, std::size_t offset) {
  if (tokens.size() < offset + 3) {
    throw std::runtime_error("Vec3 解析失败，元素数量不足");
  }
  return {std::stod(tokens[offset]), std::stod(tokens[offset + 1]),
          std::stod(tokens[offset + 2])};
}

// ============================= 配置解析器 =====================================

class ConfigParser {
 public:
  ConfigParser() = default;

  Config parse(std::istream& is) {
    Config config;
    std::unordered_map<std::string, MaterialInput> material_map;
    std::unordered_map<std::string, std::map<int, std::vector<double>>> scatter_rows;
    std::unordered_map<std::string, RegionInput> region_map;
    std::vector<std::string> region_order;

    std::string section_type;
    std::string section_name;

    std::string line;
    int line_number = 0;
    while (std::getline(is, line)) {
      ++line_number;
      auto comment_pos = line.find('#');
      if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);
      line = trim(line);
      if (line.empty()) continue;

      if (line.front() == '[' && line.back() == ']') {
        auto inside = line.substr(1, line.size() - 2);
        std::istringstream section_stream(inside);
        section_stream >> section_type;
        std::getline(section_stream, section_name);
        section_name = trim(section_name);
        if (section_type.empty()) {
          throw std::runtime_error("配置第" + std::to_string(line_number) + "行缺少节类型");
        }
        continue;
      }

      auto equal_pos = line.find('=');
      if (equal_pos == std::string::npos) {
        throw std::runtime_error("配置第" + std::to_string(line_number) + "行缺少'='");
      }
      std::string key = trim(line.substr(0, equal_pos));
      std::string value = trim(line.substr(equal_pos + 1));

      if (section_type == "global") {
        parse_global(config, key, value);
      } else if (section_type == "material") {
        auto& material = material_map[section_name];
        material.name = section_name;
        parse_material(material, scatter_rows[section_name], key, value, config.groups);
      } else if (section_type == "region") {
        auto& region = region_map[section_name];
        if (region.name.empty()) {
          region.name = section_name;
          region_order.push_back(section_name);
        }
        parse_region(region, key, value);
      } else if (section_type == "source") {
        parse_source(config.source, key, value);
        config.source.region_name = section_name;
      } else {
        throw std::runtime_error("未知配置节类型: " + section_type);
      }
    }

    // 整理材料
    for (auto& [name, mat] : material_map) {
      auto it = scatter_rows.find(name);
      if (it != scatter_rows.end()) {
        auto& rows = it->second;
        std::vector<std::vector<double>> sigma_s;
        sigma_s.resize(config.groups);
        for (std::size_t g = 0; g < config.groups; ++g) {
          auto row_it = rows.find(static_cast<int>(g));
          if (row_it == rows.end()) {
            throw std::runtime_error("材料 " + name + " 缺少 sigma_s_row" + std::to_string(g));
          }
          sigma_s[g] = row_it->second;
        }
        mat.sigma_s = std::move(sigma_s);
      }
      validate_material(mat, config.groups);
      config.materials.push_back(std::move(mat));
    }

    for (const auto& name : region_order) {
      config.regions.push_back(region_map.at(name));
    }

    if (config.groups == 0) {
      throw std::runtime_error("未在 global 节中指定 groups");
    }
    if (config.materials.empty()) {
      throw std::runtime_error("未定义任何材料");
    }
    if (config.regions.empty()) {
      throw std::runtime_error("未定义任何区域");
    }
    if (config.source.group_pdf.empty()) {
      // 若未显式指定，则默认均匀能群分布并指向外层区域
      config.source.group_pdf.assign(config.groups, 1.0 / static_cast<double>(config.groups));
      config.source.region_name = config.regions.front().name;
    }

    return config;
  }

 private:
  void parse_global(Config& config, const std::string& key, const std::string& value) {
    if (key == "groups") {
      config.groups = static_cast<std::size_t>(std::stoul(value));
    } else if (key == "histories") {
      config.options.histories = static_cast<std::size_t>(std::stoul(value));
    } else if (key == "seed") {
      config.options.seed = static_cast<unsigned>(std::stoul(value));
    } else if (key == "max_steps") {
      config.options.max_steps = static_cast<std::size_t>(std::stoul(value));
    } else {
      throw std::runtime_error("global 节中出现未知键: " + key);
    }
  }

  void parse_material(MaterialInput& material,
                      std::map<int, std::vector<double>>& scatter_rows,
                      const std::string& key, const std::string& value, std::size_t groups) {
    if (key == "sigma_a") {
      material.sigma_a = parse_double_list(value);
    } else if (key == "sigma_f") {
      material.sigma_f = parse_double_list(value);
    } else if (key == "nu") {
      material.nu = parse_double_list(value);
    } else if (key == "chi") {
      material.chi = parse_double_list(value);
    } else if (key.rfind("sigma_s_row", 0) == 0) {
      int row = std::stoi(key.substr(std::string("sigma_s_row").size()));
      scatter_rows[row] = parse_double_list(value);
    } else {
      throw std::runtime_error("material 节中出现未知键: " + key);
    }

    if (groups != 0) {
      auto check = [&](const std::vector<double>& arr, const std::string& name) {
        if (!arr.empty() && arr.size() != groups) {
          throw std::runtime_error("材料 " + material.name + " 的 " + name +
                                   " 长度与 groups 不一致");
        }
      };
      check(material.sigma_a, "sigma_a");
      check(material.sigma_f, "sigma_f");
      check(material.nu, "nu");
      check(material.chi, "chi");
    }
  }

  void parse_region(RegionInput& region, const std::string& key, const std::string& value) {
    if (key == "box") {
      auto tokens = split(value);
      if (tokens.size() != 6) {
        throw std::runtime_error("区域 " + region.name + " 的 box 需要 6 个数值");
      }
      region.min = {std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2])};
      region.max = {std::stod(tokens[3]), std::stod(tokens[4]), std::stod(tokens[5])};
    } else if (key == "material") {
      region.material_name = value;
    } else {
      throw std::runtime_error("region 节中出现未知键: " + key);
    }
  }

  void parse_source(SourceInput& source, const std::string& key, const std::string& value) {
    if (key == "group_pdf") {
      source.group_pdf = parse_double_list(value);
    } else if (key == "region") {
      source.region_name = value;
    } else {
      throw std::runtime_error("source 节中出现未知键: " + key);
    }
  }

  void validate_material(const MaterialInput& material, std::size_t groups) {
    if (material.sigma_a.size() != groups || material.sigma_f.size() != groups ||
        material.nu.size() != groups || material.chi.size() != groups) {
      throw std::runtime_error("材料 " + material.name + " 的数组维度未完成初始化");
    }
    if (material.sigma_s.size() != groups) {
      throw std::runtime_error("材料 " + material.name + " 的散射矩阵维度错误");
    }
    for (const auto& row : material.sigma_s) {
      if (row.size() != groups) {
        throw std::runtime_error("材料 " + material.name + " 的散射矩阵不是方阵");
      }
    }
  }

  // 无额外成员
};

// ============================= 求解器核心 =====================================

class RandomRaySolver {
 public:
  explicit RandomRaySolver(Config config)
      : config_(std::move(config)), rng_(config_.options.seed) {
    build_materials();
    build_geometry();
    prepare_source();
  }

  Tallies solve() {
    Tallies tallies;
    std::size_t region_count = geometry_.regions().size();
    std::size_t group_count = config_.groups;
    tallies.flux.assign(region_count, std::vector<double>(group_count, 0.0));
    tallies.absorption.assign(region_count, std::vector<double>(group_count, 0.0));
    tallies.scattering.assign(region_count, std::vector<double>(group_count, 0.0));
    tallies.fission.assign(region_count, std::vector<double>(group_count, 0.0));
    tallies.starting_histories = static_cast<double>(config_.options.histories);

    ParticleBank bank;
    for (std::size_t n = 0; n < config_.options.histories; ++n) {
      bank.push(sample_initial_particle());
    }

    const double kEps = 1e-8;
    while (!bank.empty()) {
      Particle particle = bank.pop();
      std::size_t steps = 0;

      while (particle.alive && steps < config_.options.max_steps) {
        ++steps;
        const Region* region = geometry_.locate(particle.position);
        if (!region) {
          tallies.leakage += 1.0;
          particle.alive = false;
          break;
        }

        std::size_t r_index = geometry_.region_index(region);
        const Material& mat = *region->material;
        double sigma_t = mat.sigma_t(particle.group);
        if (sigma_t <= 0.0) {
          throw std::runtime_error("材料总截面必须为正");
        }

        double dist_collision = -std::log(rng_.uniform()) / sigma_t;
        double dist_boundary = geometry_.distance_to_boundary(particle.position, particle.direction, *region);
        double flight = std::min(dist_collision, dist_boundary);

        tallies.flux[r_index][particle.group] += flight;
        tallies.absorption[r_index][particle.group] += flight * mat.sigma_a(particle.group);
        tallies.fission[r_index][particle.group] += flight * mat.sigma_f(particle.group);
        double scatter_sum = 0.0;
        for (double val : mat.scatter_row(particle.group)) scatter_sum += val;
        tallies.scattering[r_index][particle.group] += flight * scatter_sum;

        if (dist_collision < dist_boundary) {
          particle.position = particle.position + particle.direction * (flight - kEps);
          collide(mat, particle, bank);
        } else {
          particle.position = particle.position + particle.direction * (dist_boundary + kEps);
        }
      }
    }

    double normaliser = tallies.starting_histories;
    for (auto* array : {&tallies.flux, &tallies.absorption, &tallies.scattering, &tallies.fission}) {
      for (auto& region_values : *array) {
        for (double& value : region_values) value /= normaliser;
      }
    }
    tallies.leakage /= normaliser;

    return tallies;
  }

  [[nodiscard]] const Geometry& geometry() const noexcept { return geometry_; }
  [[nodiscard]] const Config& config() const noexcept { return config_; }

 private:
  void build_materials() {
    for (const auto& mat_input : config_.materials) {
      materials_.emplace(mat_input.name,
                         Material(mat_input.name, mat_input.sigma_a, mat_input.sigma_f,
                                  mat_input.nu, mat_input.sigma_s, mat_input.chi));
    }
  }

  void build_geometry() {
    geometry_ = Geometry();
    for (const auto& region_input : config_.regions) {
      auto it = materials_.find(region_input.material_name);
      if (it == materials_.end()) {
        throw std::runtime_error("区域 " + region_input.name + " 引用未知材料 " +
                                 region_input.material_name);
      }
      geometry_.add_region({region_input.name, {region_input.min, region_input.max}, &it->second});
    }
  }

  void prepare_source() {
    if (!config_.source.region_name.empty()) {
      for (const auto& region : geometry_.regions()) {
        if (region.name == config_.source.region_name) {
          source_region_index_ = geometry_.region_index(&region);
          break;
        }
      }
    }
    if (source_region_index_ == static_cast<std::size_t>(-1)) {
      source_region_index_ = 0;
    }

    if (config_.source.group_pdf.empty()) {
      source_group_pdf_.assign(config_.groups, 1.0 / static_cast<double>(config_.groups));
    } else {
      source_group_pdf_ = config_.source.group_pdf;
      double sum = std::accumulate(source_group_pdf_.begin(), source_group_pdf_.end(), 0.0);
      if (sum <= 0.0) {
        throw std::runtime_error("源项能群概率和必须为正");
      }
      for (double& value : source_group_pdf_) value /= sum;
      if (source_group_pdf_.size() != config_.groups) {
        throw std::runtime_error("源项能群概率长度与 groups 不一致");
      }
    }
  }

  Particle sample_initial_particle() {
    Particle p;
    const Region& region = geometry_.regions().at(source_region_index_);
    p.position = {region.box.min.x + rng_.uniform() * (region.box.max.x - region.box.min.x),
                  region.box.min.y + rng_.uniform() * (region.box.max.y - region.box.min.y),
                  region.box.min.z + rng_.uniform() * (region.box.max.z - region.box.min.z)};
    p.direction = normalize(rng_.isotropic_direction());
    p.group = rng_.pick(source_group_pdf_);
    return p;
  }

  void collide(const Material& mat, Particle& particle, ParticleBank& bank) {
    const auto& scatter_row = mat.scatter_row(particle.group);
    double sigma_a = mat.sigma_a(particle.group);
    double sigma_f = mat.sigma_f(particle.group);
    double scatter_sum = 0.0;
    for (double val : scatter_row) scatter_sum += val;
    double sigma_t = sigma_a + sigma_f + scatter_sum;

    double xi = rng_.uniform() * sigma_t;
    if (xi < sigma_a) {
      particle.alive = false;
      return;
    }
    xi -= sigma_a;

    if (xi < sigma_f) {
      particle.alive = false;
      spawn_fission_neutrons(mat, particle, bank);
      return;
    }
    xi -= sigma_f;

    double cumulative = 0.0;
    for (std::size_t g_out = 0; g_out < scatter_row.size(); ++g_out) {
      cumulative += scatter_row[g_out];
      if (xi <= cumulative) {
        particle.group = g_out;
        break;
      }
    }
    particle.direction = normalize(rng_.isotropic_direction());
  }

  void spawn_fission_neutrons(const Material& mat, const Particle& parent, ParticleBank& bank) {
    double nu = mat.nu(parent.group);
    int integer = static_cast<int>(std::floor(nu));
    double fractional = nu - integer;
    int offspring = integer;
    if (rng_.uniform() < fractional) ++offspring;
    if (offspring == 0) return;

    const auto& chi = mat.chi();
    for (int i = 0; i < offspring; ++i) {
      Particle child;
      child.position = parent.position;
      child.direction = normalize(rng_.isotropic_direction());
      child.group = rng_.pick(chi);
      bank.push(child);
    }
  }

 private:
  Config config_;
  RandomEngine rng_;
  std::unordered_map<std::string, Material> materials_;
  Geometry geometry_;
  std::size_t source_region_index_ = static_cast<std::size_t>(-1);
  std::vector<double> source_group_pdf_;
};

}  // namespace random_ray_solver

// ============================= 主程序入口 =====================================

int main(int argc, char** argv) {
  using namespace random_ray_solver;

  std::string input_path = "config/example_reactor.rr";
  std::size_t override_histories = 0;
  unsigned override_seed = 0;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--input" && i + 1 < argc) {
      input_path = argv[++i];
    } else if (arg == "--histories" && i + 1 < argc) {
      override_histories = static_cast<std::size_t>(std::stoul(argv[++i]));
    } else if (arg == "--seed" && i + 1 < argc) {
      override_seed = static_cast<unsigned>(std::stoul(argv[++i]));
    } else if (arg == "--help") {
      std::cout << "Usage: ./random_ray_solver [--input path] [--histories N] [--seed S]\n";
      return 0;
    } else {
      std::cerr << "未知命令行参数: " << arg << "\n";
      return 1;
    }
  }

  std::ifstream ifs(input_path);
  if (!ifs) {
    std::cerr << "无法打开输入文件: " << input_path << "\n";
    return 1;
  }

  ConfigParser parser;
  Config config;
  try {
    config = parser.parse(ifs);
  } catch (const std::exception& e) {
    std::cerr << "配置解析错误: " << e.what() << "\n";
    return 1;
  }

  if (override_histories > 0) config.options.histories = override_histories;
  if (override_seed > 0) config.options.seed = override_seed;

  try {
    RandomRaySolver solver(std::move(config));
    Tallies results = solver.solve();
    const auto& regions = solver.geometry().regions();

    std::cout << "================ 随机射线多群输运统计 ================\n";
    std::cout << "初始历史数: " << results.starting_histories << "\n";
    std::cout << "几何区域数: " << regions.size() << "\n";
    std::cout << "能群数: " << solver.config().groups << "\n";
    std::cout << "漏出归一化权重: " << std::setprecision(6) << results.leakage << "\n";
    std::cout << "-------------------------------------------------------\n";

    for (std::size_t r = 0; r < regions.size(); ++r) {
      std::cout << "区域: " << regions[r].name << " (材料: " << regions[r].material->name() << ")\n";
      for (std::size_t g = 0; g < results.flux[r].size(); ++g) {
        std::cout << "  能群 " << g << " 通量=" << std::setw(12) << results.flux[r][g]
                  << " 吸收率=" << std::setw(12) << results.absorption[r][g]
                  << " 散射率=" << std::setw(12) << results.scattering[r][g]
                  << " 裂变率=" << std::setw(12) << results.fission[r][g] << '\n';
      }
      std::cout << "-------------------------------------------------------\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "求解器运行失败: " << e.what() << "\n";
    return 1;
  }

  return 0;
}

