#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// TrackLengthTally 展示了“策略式 Tallier + 自动单元归档”的思想：
// - Segment 数据结构清晰表达粒子路径。
// - TallyPolicy 抽象 tallies 需求，可扩展到能群/角度。
// - TrackLengthTally 聚合 policy 输出，展示组合与缓存。

struct Segment {
  std::string cell_name;
  double path_length_cm;
  double flux_weight;
};

class TallyPolicy {
public:
  virtual ~TallyPolicy() = default;
  virtual std::string bin_label(const Segment& seg) const = 0;
  virtual double contribution(const Segment& seg) const = 0;
};

class FluxPolicy final : public TallyPolicy {
public:
  std::string bin_label(const Segment& seg) const override { return seg.cell_name; }

  double contribution(const Segment& seg) const override
  {
    return seg.path_length_cm * seg.flux_weight;
  }
};

class TrackLengthTally {
public:
  explicit TrackLengthTally(std::unique_ptr<TallyPolicy> policy)
    : policy_{std::move(policy)}
  {
  }

  void accumulate(const Segment& seg)
  {
    const std::string key = policy_->bin_label(seg);
    tallies_[key] += policy_->contribution(seg);
  }

  void print() const
  {
    std::cout << std::fixed << std::setprecision(4);
    for (const auto& [label, value] : tallies_) {
      std::cout << label << ": " << value << '\n';
    }
  }

private:
  std::unique_ptr<TallyPolicy> policy_;
  std::map<std::string, double> tallies_;
};

int main()
{
  TrackLengthTally tally{std::make_unique<FluxPolicy>()};

  std::vector<Segment> history{
    {"fuel", 0.5, 1.0},
    {"moderator", 1.5, 0.9},
    {"fuel", 0.3, 1.1},
  };

  for (const auto& seg : history) {
    tally.accumulate(seg);
  }

  tally.print();
}
