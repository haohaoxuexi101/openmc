#include <iostream>
#include <vector>

struct Particle {
  double energy;
  double weight;
};

class SecondaryBank {
public:
  void add(const Particle& p) { bank_.push_back(p); }

  void split(Particle& parent, int n)
  {
    double new_weight = parent.weight / (n + 1);
    parent.weight = new_weight;
    for (int i = 0; i < n; ++i) {
      bank_.push_back(Particle{parent.energy * 0.5, new_weight});
    }
  }

  void dump() const
  {
    std::cout << "Secondary bank size=" << bank_.size() << '\n';
    for (const auto& p : bank_) {
      std::cout << "  energy=" << p.energy << " weight=" << p.weight << '\n';
    }
  }

private:
  std::vector<Particle> bank_;
};

int main()
{
  Particle parent{2.0e5, 1.0};
  SecondaryBank bank;
  bank.split(parent, 3);
  bank.add(Particle{1.0e5, 0.2});
  std::cout << "Parent new weight=" << parent.weight << '\n';
  bank.dump();
}
