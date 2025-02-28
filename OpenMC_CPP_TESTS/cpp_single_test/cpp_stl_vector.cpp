#include <iostream>
#include <vector>

int main()
{
  std::vector<int> test_vector {4, 8, 5, 6};
  test_vector.push_back(2);
  test_vector.pop_back();
  test_vector.insert(test_vector.begin() + 2, 3);
  // basic print method
  for (int i = 0; i < test_vector.size(); i++) {
    std::cout << test_vector[i] << " ";
  }
  std::cout << std::endl;

  test_vector.clear();
  test_vector.reserve(1);
  test_vector.push_back(1);
  test_vector.push_back(2);
  test_vector.push_back(3);
  test_vector.push_back(4);
  test_vector.erase(test_vector.begin() + 1);

  std::cout << "the end: " << test_vector.back() << std::endl;
  std::cout << "the begin: " << test_vector.front() << std::endl;

  // another print method
  for (const auto& element : test_vector) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
  return 0;
}