#include "vector.hpp"
#include <cmath>
#include <iostream>

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
vector operator*(vector a, vector b) {
  return vector(a.x * b.x, a.y * b.y, a.z * b.y);
}
vectorf vecToVecf(vector a) { return vectorf(a.x, a.y, a.z); }
vector vecfToVec(vectorf a) { return vector(a.x, a.y, a.z); }
double vector::mag() { return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)); }

vector vector::unit() { return *this / mag(); }
vector vector::opposite() { return vector(-x, -y, -z); }
vector vector::transform_return(vector ihat_prime, vector jhat_prime) {
  return vector(ihat_prime.x * x + jhat_prime.x * y,
                ihat_prime.y * x + jhat_prime.y * y);
}
void vector::transform(vector ihat_prime, vector jhat_prime) {
  double setx = ihat_prime.x * x + jhat_prime.x * y;
  y = ihat_prime.y * x + jhat_prime.y * y;
  x = setx;
}

void printvec(vector a) {
  std::cout << "vector: {" << a.x << ", " << a.y << ", " << a.z << "}\n";
}
void printvec(vector a, std::string name) {
  std::cout << name << ": {" << a.x << ", " << a.y << ", " << a.z << "}\n";
}
