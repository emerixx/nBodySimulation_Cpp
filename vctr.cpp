#include "vctr.h"
#include <cmath>
#include <iostream>
#include <string>

vector operator-(vector a, vector b) {
  return vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

vector operator+(vector a, vector b) {
  return vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

vector operator/(vector a, double b) {
  return vector(a.x / b, a.y / b, a.z / b);
}

vector operator*(vector a, double b) {
  return vector(a.x * b, a.y * b, a.z * b);
}

vectorf vecToVecf(vector a) { return vectorf(a.x, a.y, a.z); }

double magnitude(vector a) {
  return sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2));
}

vector normalize(vector a) {
  double mag = magnitude(a);
  return vector(a.x / mag, a.y / mag, a.z / mag);
}
vector opposite(vector a) { return vector(-a.x, -a.y, -a.z); }
void print(vector a) {
  std::cout << "vec.x: " + std::to_string(a.x) << "\n";
  std::cout << "vec.y: " + std::to_string(a.y) << "\n";
  std::cout << "vec.z: " + std::to_string(a.z) << "\n";
}
