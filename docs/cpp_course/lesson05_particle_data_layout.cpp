#include <cmath>
#include <iostream>
#include <string>

// 通过 ParticleData/Particle 的分层封装展示 OpenMC 如何既保持内存紧凑又提供面向对象接口。

struct ParticleData {
  double position[3];
  double direction[3];
  double speed;
  double energy;
};

class Particle {
public:
  explicit Particle(ParticleData& data, std::string label)
    : data_{data}
    , label_{std::move(label)}
  {
  }

  void move(double distance)
  {
    for (int i = 0; i < 3; ++i) {
      data_.position[i] += distance * data_.direction[i];
    }
  }

  void scatter(double mu)
  {
    data_.direction[2] = mu;
    const double perp = std::sqrt(1.0 - mu * mu);
    data_.direction[0] = perp;
    data_.direction[1] = 0.0;
  }

  void print_state() const
  {
    std::cout << label_ << " pos=(" << data_.position[0] << ',' << data_.position[1] << ','
              << data_.position[2] << ") energy=" << data_.energy << '\n';
  }

private:
  ParticleData& data_;
  std::string label_;
};

int main()
{
  ParticleData data{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.5, 1.0e5};
  Particle particle{data, "primary"};
  particle.print_state();
  particle.move(3.0);
  particle.scatter(0.6);
  particle.move(1.0);
  particle.print_state();
}
