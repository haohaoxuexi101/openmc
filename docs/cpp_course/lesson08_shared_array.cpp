#include <atomic>
#include <iostream>
#include <thread>
#include <vector>

class SharedArray {
public:
  explicit SharedArray(std::size_t capacity) : data_(capacity) {}

  void push(double value)
  {
    auto index = write_index_.fetch_add(1, std::memory_order_relaxed);
    if (index < data_.size()) {
      data_[index] = value;
    }
  }

  std::vector<double> snapshot() const
  {
    return std::vector<double>(data_.begin(), data_.begin() + write_index_.load());
  }

private:
  std::vector<double> data_;
  std::atomic<std::size_t> write_index_{0};
};

int main()
{
  SharedArray bank(16);
  auto worker = [&bank](double base) {
    for (int i = 0; i < 4; ++i) {
      bank.push(base + i);
    }
  };

  std::thread t1(worker, 0.0);
  std::thread t2(worker, 100.0);
  std::thread t3(worker, 200.0);

  t1.join();
  t2.join();
  t3.join();

  auto values = bank.snapshot();
  std::cout << "Collected " << values.size() << " values\n";
  for (double v : values) {
    std::cout << "  " << v << '\n';
  }
}
