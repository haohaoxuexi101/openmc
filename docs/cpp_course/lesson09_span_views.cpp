#include <cstddef>
#include <iostream>
#include <numeric>
#include <vector>

template <typename T>
class Span {
public:
  Span(T* data, std::size_t size) : data_{data}, size_{size} {}

  T& operator[](std::size_t idx) { return data_[idx]; }
  const T& operator[](std::size_t idx) const { return data_[idx]; }

  T* begin() { return data_; }
  T* end() { return data_ + size_; }
  const T* begin() const { return data_; }
  const T* end() const { return data_ + size_; }

  std::size_t size() const { return size_; }

private:
  T* data_;
  std::size_t size_;
};

int main()
{
  std::vector<double> xs{1.0, 2.0, 3.0, 4.0};
  Span<double> tail(xs.data() + 1, xs.size() - 1);
  for (double& value : tail) {
    value *= 2.0;
  }

  double sum = std::accumulate(xs.begin(), xs.end(), 0.0);
  std::cout << "Sum with span-modified tail = " << sum << '\n';
}
