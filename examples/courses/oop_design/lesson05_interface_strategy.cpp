#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

// 课程 5：策略模式与接口稳定性
// 演示如何使用抽象基类定义计分策略，便于扩展新的响应函数，
// 对应 OpenMC 中 TallyFilter / Score 的策略化设计。

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
