#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <string>
#include <vector>
using vector_container = std::vector<std::vector<double>>;
using vector_double = std::vector<double>;
namespace mml {
class matrix;
class vector2 {
public:
  double x;
  double y;
  vector2(double x_in, double y_in) : x(x_in), y(y_in) {}
  vector2() : x(0), y(0) {}
  vector2 unit();
  double mag();
  vector2 opposite();
  void transform(vector2 ihat_prime, vector2 jhat_prime);
  void transform(matrix &transformation_matrix);
  vector2 transform_return(vector2 ihat_prime, vector2 jhat_prime);
  vector2 transform_return(matrix &transformation_matrix);
  vector2 abs();
};

double dot(vector2 a, vector2 b);

void printvec(vector2 a);
void printvec(vector2 a, std::string name);

class matrix {
private:
  uint max_size = 32; // defined to prevent accidents
  vector_container matrx_vec;

  mml::vector2 matrix_size;

public:
  matrix(mml::vector2 arg_size);

  void print();
  mml::vector2 getSize();
  double determinant();
  void invert(); // turn this matrix to its own inverse, A = A^-1
  mml::matrix inverse();
  double &operator[](size_t row, size_t col);
  vector_double &operator[](size_t row);
};

} // namespace mml

mml::vector2 operator-(mml::vector2 a, mml::vector2 b);
mml::vector2 operator+(mml::vector2 a, mml::vector2 b);
mml::vector2 operator/(mml::vector2 a, double b);
mml::vector2 operator*(mml::vector2 a, double b);
mml::vector2 operator*(mml::vector2 a, mml::vector2 b);
mml::vector2 operator/(mml::vector2 a, mml::vector2 b);
bool operator>(mml::vector2 a, mml::vector2 b);
bool operator>=(mml::vector2 a, mml::vector2 b);
bool operator<(mml::vector2 a, mml::vector2 b);
bool operator<=(mml::vector2 a, mml::vector2 b);
bool operator==(mml::vector2 a, mml::vector2 b);
bool operator!=(mml::vector2 a, mml::vector2 b);
/*
class vector3 {
public:
  double x;
  double y;
  double z;
  vector3(double x_in, double y_in) : x(x_in), y(y_in), z(0) {}
  vector3(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in) {}
  vector3() : x(0), y(0), z(0) {}
  vector3 unit();
  double mag();
  vector3 opposite();
  vector3 abs();
  vector3 set_z_temp(double a);
};

struct vectorf {
  float x;
  float y;
  float z;
  vectorf(float x_in, float y_in) : x(x_in), y(y_in), z(0) {}
  vectorf(float x_in, float y_in, float z_in) : x(x_in), y(y_in), z(z_in) {}
  vectorf() : x(0), y(0), z(0) {}
};
*/
double sign(double a);

#endif // VECTOR_H
