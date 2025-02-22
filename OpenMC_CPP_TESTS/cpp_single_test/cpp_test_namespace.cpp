#include <iostream>
using namespace std;
void swap(int& a, int& b)
{
  int temp = a;
  a = b;
  b = temp;
}
int main()
{
  int a = 1;
  int b = 2;
  swap(a, b);
  cout << "a = " << a << ", b = " << b << "\n";
  cout << "Hello, World!" << "\n";
  return 0;
}