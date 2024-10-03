#include "vector.hpp"
#include <cmath>
#include <iostream>
mml::vector2 operator-(mml::vector2 a, mml::vector2 b) {
  return mml::vector2(a.x - b.x, a.y - b.y);
}
mml::vector2 operator+(mml::vector2 a, mml::vector2 b) {
  return mml::vector2(a.x + b.x, a.y + b.y);
}

mml::vector2 operator/(mml::vector2 a, double b) {
  return mml::vector2(a.x / b, a.y / b);
}
mml::vector2 operator*(mml::vector2 a, double b) {
  return mml::vector2(a.x * b, a.y * b);
}
mml::vector2 operator*(mml::vector2 a, mml::vector2 b) {
  return mml::vector2(a.x * b.x, a.y * b.y);
}
mml::vector2 operator/(mml::vector2 a, mml::vector2 b) {
  return mml::vector2(a.x / b.x, a.y / b.y);
}
bool operator>(mml::vector2 a, mml::vector2 b) {
  return a.x > b.x && a.y > b.y;
}
bool operator>=(mml::vector2 a, mml::vector2 b) {
  return a.x >= b.x && a.y >= b.y;
}
bool operator<(mml::vector2 a, mml::vector2 b) {
  return a.x < b.x && a.y < b.y;
}
bool operator<=(mml::vector2 a, mml::vector2 b) {
  return a.x <= b.x && a.y <= b.y;
}
bool operator==(mml::vector2 a, mml::vector2 b) {
  return a.x == b.x && a.y == b.y;
}
bool operator!=(mml::vector2 a, mml::vector2 b) {
  return a.x != b.x || a.y != b.y;
}

double mml::dot(mml::vector2 a, mml::vector2 b) {
  // source:
  // https://math.libretexts.org/Bookshelves/Calculus/Calculus_3e_(Apex)/10%3A_Vectors/10.03%3A_The_Dot_Product
  return a.x * b.x + a.y * b.y;
}

double mml::vector2::mag() { return sqrt(pow(x, 2) + pow(y, 2)); }

mml::vector2 mml::vector2::unit() { return *this / mag(); }
mml::vector2 mml::vector2::opposite() { return mml::vector2(-x, -y); }
mml::vector2 mml::vector2::transform_return(mml::vector2 ihat_prime,
                                            mml::vector2 jhat_prime) {
  return mml::vector2(ihat_prime.x * x + jhat_prime.x * y,
                      ihat_prime.y * x + jhat_prime.y * y);
}
mml::vector2
mml::vector2::transform_return(mml::matrix &transformation_matrix) {
  if (transformation_matrix.getSize() != mml::vector2(2, 2)) {
    std::cout << "umm, matrix is too big for me UwU (error transforming a "
                 "vector, matrix is bigger than 2x2)\n";
    return vector2(0, 0);
  }
  mml::matrix mtrx = transformation_matrix;
  mml::vector2 ihat_prime = {mtrx[0, 0], mtrx[0, 1]};
  mml::vector2 jhat_prime = {mtrx[1, 0], mtrx[1, 1]};
  return mml::vector2(ihat_prime.x * x + jhat_prime.x * y,
                      ihat_prime.y * x + jhat_prime.y * y);
}
void mml::vector2::transform(mml::vector2 ihat_prime, mml::vector2 jhat_prime) {
  double setx = ihat_prime.x * x + jhat_prime.x * y;
  y = ihat_prime.y * x + jhat_prime.y * y;
  x = setx;
}
void mml::vector2::transform(mml::matrix &transformation_matrix) {
  if (transformation_matrix.getSize() != mml::vector2(2, 2)) {
    std::cout << "umm, matrix is too big for me UwU (error transforming a "
                 "vector, matrix is bigger than 2x2)\n";
    return;
  }
  mml::matrix mtrx = transformation_matrix;
  mml::vector2 ihat_prime = {mtrx[0, 0], mtrx[1, 0]};
  mml::vector2 jhat_prime = {mtrx[0, 1], mtrx[1, 1]};
  double setx = ihat_prime.x * x + jhat_prime.x * y;
  y = ihat_prime.y * x + jhat_prime.y * y;
  x = setx;
}
mml::vector2 mml::vector2::abs() {
  return mml::vector2(std::abs(x), std::abs(y));
};
void mml::printvec(mml::vector2 a) {
  std::cout << "vector: {" << a.x << ", " << a.y << "}\n";
}
void mml::printvec(mml::vector2 a, std::string name) {
  std::cout << name << ": {" << a.x << ", " << a.y << "}\n";
}
// MATRIX --------------
mml::matrix::matrix(mml::vector2 arg_size) {
  if (arg_size.x > max_size || arg_size.y > max_size) {
    std::cout
        << "ERROR: trying to initialize an matrix thats bigger than max_size, "
           "edit matrix.hpp if this wasnt a mistake\nexiting...\n";
  }
  matrix_size = arg_size.transform_return({0, 1}, {1, 0});
  for (int row = 0; row < matrix_size.y; row++) {
    matrx_vec.push_back({0});
    for (int col = 0; col < matrix_size.x; col++) {
      matrx_vec[row].push_back(0);
      matrx_vec[row][col] = 0;
    }
  }
}
void mml::matrix::print() {
  std::cout << "----\n";
  for (int row = 0; row < matrix_size.y; row++) {
    if (row > 0) {
      std::cout << "\n";
    }
    std::cout << "|";
    for (int col = 0; col < matrix_size.x; col++) {
      std::cout << matrx_vec[row][col];
      if (col < matrix_size.x - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "|";
  }
  std::cout << "\n----\n";
}
mml::vector2 mml::matrix::getSize() { return matrix_size; }
double &mml::matrix::operator[](size_t row, size_t col) {

  if (row >= matrix_size.y || col >= matrix_size.x) {
    std::cout
        << "ERROR: trying to access element outside of matrix\nexiting...\n";
    exit(-1);
  }
  return matrx_vec[row][col];
}

vector_double &mml::matrix::operator[](size_t row) {
  if (row >= matrix_size.y) {
    std::cout
        << "ERROR: trying to access element outside of matrix\nexiting...\n";
    exit(-1);
  }
  return matrx_vec[row];
}

double mml::matrix::determinant() {
  if (getSize() != vector2(2, 2)) {
    std::cout << "umm, matrix is too big for me UwU (error finding the "
                 "determinant of a "
                 "matrix, matrix is bigger than 2x2)\n";
    return 0;
  }
  matrix A = *this;
  // det(A) = ad-bc
  // det(A) = [0,0]*[1,1] - [1,0]*[0,1]
  return A[0, 0] * A[1, 1] - A[1, 0] * A[0, 1];
}
void mml::matrix::invert() {
  if (getSize() != vector2(2, 2)) {
    std::cout << "umm, matrix is too big for me UwU (error finding the "
                 "inverse matrix of a "
                 "matrix, matrix is bigger than 2x2)\n";
    return;
  }
  matrix A = *this;
  matrix temp_mat = mml::matrix({2, 2});
  double det = A.determinant();
  // B:
  //  --         --
  //  |  A.d -A.b |
  //  | -A.c  A.a |
  //  --         --
  // A^-1 = det(A)^-1 * B
  temp_mat[0, 0] = 1 / det * A[1, 1];
  temp_mat[0, 1] = 1 / det * -A[0, 1];
  temp_mat[1, 0] = 1 / det * -A[1, 0];
  temp_mat[1, 1] = 1 / det * A[0, 0];
  *this = temp_mat;
}
mml::matrix mml::matrix::inverse() {
  if (getSize() != vector2(2, 2)) {
    std::cout << "umm, matrix is too big for me UwU (error finding the "
                 "inverse matrix of a "
                 "matrix, matrix is bigger than 2x2)\n";
    return mml::matrix({2, 2});
  }
  matrix A = *this;
  matrix temp_mat = mml::matrix({2, 2});
  double det = A.determinant();
  // B:
  //  --         --
  //  |  A.d -A.b |
  //  | -A.c  A.a |
  //  --         --
  // A^-1 = det(A)^-1 * B
  temp_mat[0, 0] = 1 / det * A[1, 1];
  temp_mat[0, 1] = 1 / det * -A[0, 1];
  temp_mat[1, 0] = 1 / det * -A[1, 0];
  temp_mat[1, 1] = 1 / det * A[0, 0];
  return temp_mat;
}

double sign(double a) { return a / std::abs(a); }
