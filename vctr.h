#ifndef VCTR_H
#define VCTR_H

struct vector {
  double x;
  double y;
  double z;
  vector(double x_in, double y_in) : x(x_in), y(y_in), z(0) {}
  vector(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in) {}
  vector() : x(0), y(0), z(0) {}
};

struct vectorf {
  float x;
  float y;
  float z;
  vectorf(float x_in, float y_in) : x(x_in), y(y_in), z(0) {}
  vectorf(float x_in, float y_in, float z_in) : x(x_in), y(y_in), z(z_in) {}
  vectorf() : x(0), y(0), z(0) {}
};

vector operator-(vector a, vector b);
vector operator+(vector a, vector b);
vector operator/(vector a, double b);
vector operator*(vector a, double b);
vectorf vecToVecf(vector a);
double magnitude(vector a);
vector normalize(vector a);
vector opposite(vector a);
void print(vector a);

#endif // VCTR_H
