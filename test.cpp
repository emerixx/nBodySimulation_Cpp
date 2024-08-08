#include <iostream>

struct sa {
  int x;
  int y;
  int z;
};
struct sb {
  int w;
  int k;
  int j;
};
int operator+(sa a, sb b) { return 0; }
int operator+(sb a, sa b) { return 1; }
int main() {

  sa m = {1, 2, 3};
  sb n = {4, 5, 6};
  std::cout << m + n << "\n";
  std::cout << n + m << "\n";
  std::cout << "hello world\n";
  return 0;
}
