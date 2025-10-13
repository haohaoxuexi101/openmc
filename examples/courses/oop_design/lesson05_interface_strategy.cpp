#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

// 课程 5：策略模式与接口稳定性
// 借鉴说明：该模式来自 OpenMC `src/tallies/tally.h` 中 Score 组合 TallyFilter
// 的结构——通过抽象接口定义可扩展的计分行为。示例保留“策略对象 + Tally”
// 的协作方式，但将计分函数压缩为简单数学模型，用以强调接口稳定性。

class ScoreStrategy {
public:
  virtual ~ScoreStrategy() = default;
  virtual double score(double flux, double xs) const = 0;
  virtual std::string description() const = 0;
};

class ReactionRateScore : public ScoreStrategy {
public:
  double score(double flux, double xs) const override { return flux * xs; }
  std::string description() const override { return "反应率计分"; }
};

class AbsorptionProbabilityScore : public ScoreStrategy {
public:
  double score(double flux, double xs) const override { return 1.0 - std::exp(-flux * xs); }
  std::string description() const override { return "吸收概率计分"; }
};

class Tally {
public:
  explicit Tally(std::shared_ptr<ScoreStrategy> strategy)
      : strategy_(std::move(strategy)) {}

  void evaluate(double flux, double xs) const {
    std::cout << strategy_->description() << " -> 结果: "
              << strategy_->score(flux, xs) << "\n";
  }

private:
  std::shared_ptr<ScoreStrategy> strategy_;
};

int main() {
  auto reaction_rate = std::make_shared<ReactionRateScore>();
  auto absorption = std::make_shared<AbsorptionProbabilityScore>();

  Tally tally_rr(reaction_rate);
  Tally tally_abs(absorption);

  double flux = 0.8;  // 单位化通量
  double sigma = 0.5; // 有效截面

  tally_rr.evaluate(flux, sigma);
  tally_abs.evaluate(flux, sigma);
  return 0;
}
