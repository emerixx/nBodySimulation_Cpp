#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <raylib.h>
#include <rlgl.h>
#include <string>

const double pi = 3.1415926535897932384626433832795028841971693993751;
const double G = 4 * pow(pi, 2); // AU^3 MO^-1 Year^-2
const double constTimeStep = pow(10, -4);
const int nOfBodies = 2;
const bool hasDE = false;
double speed = 100;
double B[7][6];

double CH[] = {0,         16.0 / 135, 0, 6656.0 / 12825, 28561.0 / 56430,
               -9.0 / 50, 2.0 / 55};
void print(std::string out) { std::cout << out << "\n"; }
std::string str(double num) { return std::to_string(num); }

struct vector {
  double x;
  double y;
  double z;
  vector(double x_in, double y_in) : x(x_in), y(y_in), z(0) {}
  vector(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in) {}
  vector() : x(0), y(0), z(0) {}
  double magnitude() { return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)); }
  vector unit() {
    double mag = magnitude();
    return vector(x / mag, y / mag, z / mag);
  }
  vector opposite() { return vector(-x, -y, -z); }
  void Print() {
    print("vector.x: " + str(x));
    print("vector.y: " + str(y));
    print("vector.z: " + str(z));
  }
};

struct vectorf {
  float x;
  float y;
  float z;
  vectorf(float x_in, float y_in) : x(x_in), y(y_in), z(0) {}
  vectorf(float x_in, float y_in, float z_in) : x(x_in), y(y_in), z(z_in) {}
  vectorf() : x(0), y(0), z(0) {}
};
struct bs_struct {
  vector position;
  vector velocity;
  double mass;

  bs_struct(vector p, vector v, double m) : position(p), velocity(v), mass(m) {}
  bs_struct() : position({0, 0}), velocity({0, 0}), mass(0){};
  void Print() {
    print("PRINT bs_struct: START--------------------");
    print("position:");
    position.Print();
    print("velocity:");
    velocity.Print();
    print("mass: " + str(mass));
    print("PRINT bs_struct: END----------------------");
  }
};

struct ds_struct {
  vector velocity;
  vector acceleration;

  ds_struct(vector v, vector a) : velocity(v), acceleration(a) {}
  ds_struct() : velocity({0, 0}), acceleration({0, 0}){};

  void Print() {
    print("PRINT ds_struct: START--------------------");
    print("velocity:");
    velocity.Print();
    print("acceleration:");
    acceleration.Print();
    print("PRINT ds_struct: END----------------------");
  }
};

Color bodyColor[] = {RED, BLUE, GREEN};
vectorf winSize = vectorf(1600, 900);
vectorf winCenter = vectorf(winSize.x / 2, winSize.y / 2);
double drawScale = 50;

int circleSegments = 32;

vectorf camera_startPos = vectorf{0, 0, 250};
vectorf camera_startTarget = vectorf(0, 0, 0);
vectorf camera_up = vectorf(0, 1, 0);
float camera_startFovy = 60;

Camera camera;
Image img;
Texture texture;
////////////////////////
//   Setup 2 bodies   //
////////////////////////

std::array<bs_struct, nOfBodies> g2b(std::array<bs_struct, nOfBodies> cbs) {
  cbs[0] = bs_struct({-1, 0}, {0, -1}, 1);
  cbs[1] = bs_struct({1, 0}, {0, 1}, 1);
  return cbs;
}
std::array<bs_struct, nOfBodies> g3b(std::array<bs_struct, nOfBodies> cbs) {
  cbs[0] = bs_struct({-1, 0}, {0, -1}, 1);
  cbs[1] = bs_struct({1, 0}, {0, 1}, 1);
  cbs[2] = bs_struct({0, 1}, {0, 0}, 1);
  return cbs;
}

void drawGrid(int slices, float spacing, Color color) {

  int halfSlices = slices / 2;

  rlBegin(RL_LINES);
  for (int i = -halfSlices; i <= halfSlices; i++) {

    rlColor3f(color.r / 255.f, color.g / 255.f, color.b / 255.f);

    rlVertex3f((float)i * spacing, 0.0f, (float)-halfSlices * spacing);
    rlVertex3f((float)i * spacing, 0.0f, (float)halfSlices * spacing);

    rlVertex3f((float)-halfSlices * spacing, 0.0f, (float)i * spacing);
    rlVertex3f((float)halfSlices * spacing, 0.0f, (float)i * spacing);
  }
  rlEnd();
}

void drawCircleOld(vector centerOffset, double radius, Color color) {
  int segments = circleSegments;
  double step = 2 * pi / segments;
  double angle = 0;
  rlBegin(RL_TRIANGLES);
  rlSetTexture(texture.id);
  for (int i = 0; i < segments; i++) {
    rlColor4ub(color.r, color.g, color.b, color.a);
    rlVertex2f(0, 0);
    rlVertex2f(sin(angle + step) * radius, cos(angle + step) * radius);
    rlVertex2f(sin(angle) * radius, cos(angle) * radius);
    angle += step;
  }
  rlEnd();
}
void drawCircle(vectorf centerOffset, float radius, Color color, float scale) {
  int segments = circleSegments;
  double step = 2 * pi / segments;
  double angle = 0;
  vectorf center = {centerOffset.x * scale, centerOffset.y * scale};
  vectorf vertex2 = {sinf(angle + step) * (float)radius,
                     cosf(angle + step) * (float)radius};
  vectorf vertex2old = {0, (float)radius};
  rlBegin(RL_TRIANGLES);
  rlSetTexture(texture.id);
  for (int i = 0; i < segments; i++) {

    rlColor4ub(color.r, color.g, color.b, color.a);
    rlVertex2f(center.x, center.y);
    rlVertex2f(vertex2.x + center.x, vertex2.y + center.y);
    rlVertex2f(vertex2old.x + center.x, vertex2old.y + center.y);
    angle += step;
    vertex2old = vertex2;
    vertex2 = {sinf(angle + step) * (float)radius,
               cosf(angle + step) * (float)radius};
  }
  rlEnd();
}

vector subtractVec(vector vec1, vector vec2) {
  return vector(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
}
vector addVec(vector v1, vector v2) {
  return vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}
vector divideVec(vector v, double s) {
  return vector(v.x / s, v.y / s, v.z / s);
}
vector multiplyVec(vector v, double s) {
  return vector(v.x * s, v.y * s, v.z * s);
}

vectorf vecToVecf(vector v) { return vectorf(v.x, v.y); }

std::array<ds_struct, nOfBodies>
deltaStateFunction(std::array<bs_struct, nOfBodies> bsLoc) {
  vector ForceArr[nOfBodies][nOfBodies];
  vector Force[nOfBodies];
  for (int i = 0; i < nOfBodies; i++) {
    for (int j = 0; j < nOfBodies; j++) {
      if (i == j)
        continue;
      if (i > j) {
        ForceArr[i][j] = ForceArr[j][i].opposite();
        continue;
      }
      vector distance = subtractVec(bsLoc[j].position, bsLoc[i].position);
      double distanceMagnitude = distance.magnitude();
      vector direction = distance.unit();
      double forceMagnitude = 0;
      if (distanceMagnitude > 0) {
        forceMagnitude =
            (G * bsLoc[i].mass * bsLoc[j].mass) / pow(distanceMagnitude, 2);
      }
      ForceArr[i][j] = multiplyVec(direction, forceMagnitude);
    }
    for (int j = 0; j < nOfBodies; j++) {
      if (i == j)
        continue;
      Force[i] = addVec(Force[i], ForceArr[i][j]);
    }
  }
  std::array<ds_struct, nOfBodies> out;
  for (int i = 0; i < nOfBodies; i++) {
    out[i] = ds_struct(bsLoc[i].velocity, divideVec(Force[i], bsLoc[i].mass));
  }
  return out;
}

std::array<bs_struct, nOfBodies>
multiplyDsByTime(std::array<ds_struct, nOfBodies> ds, double t) {
  std::array<bs_struct, nOfBodies> out;
  for (int i = 0; i < nOfBodies; i++) {
    out[i] = bs_struct(multiplyVec(ds[i].velocity, t),
                       multiplyVec(ds[i].acceleration, t), 0);
  }
  return out;
}
std::array<bs_struct, nOfBodies>
multiplyBsByScalar(std::array<bs_struct, nOfBodies> bs, double scalar) {
  std::array<bs_struct, nOfBodies> out;
  for (int i = 0; i < nOfBodies; i++) {
    out[i] =
        bs_struct(multiplyVec(bs[i].position, scalar),
                  multiplyVec(bs[i].velocity, scalar), bs[i].mass * scalar);
  }
  return out;
}
std::array<bs_struct, nOfBodies> addBs(std::array<bs_struct, nOfBodies> bs1,
                                       std::array<bs_struct, nOfBodies> bs2) {
  std::array<bs_struct, nOfBodies> out;
  for (int i = 0; i < nOfBodies; i++) {
    out[i] = bs_struct(addVec(bs1[i].position, bs2[i].position),
                       addVec(bs1[i].velocity, bs2[i].velocity),
                       bs1[i].mass + bs2[i].mass);
  }
  return out;
}
std::array<bs_struct, nOfBodies> RKF45(std::array<bs_struct, nOfBodies> cbs) {
  double dt = constTimeStep * GetFrameTime() * speed;
  if (!hasDE) {
    dt = constTimeStep;
  }
  std::array<std::array<bs_struct, nOfBodies>, 7> k;

  for (int i = 1; i < 7; i++) {
    std::array<bs_struct, nOfBodies> dsfArg = cbs;
    for (int j = 1; j < i; j++) {
      dsfArg = addBs(dsfArg, multiplyBsByScalar(k[j], B[i][j]));
    }
    k[i] = multiplyDsByTime(deltaStateFunction(dsfArg), dt);
  }
  /*
  k[1] = multiplyDsByTime(deltaStateFunction(cbs), dt);
  k[2] = multiplyDsByTime(
      deltaStateFunction(addBs(cbs, multiplyBsByScalar(k[1], B[2][1]))), dt);
  k[3] = multiplyDsByTime(
      deltaStateFunction(addBs(addBs(cbs, multiplyBsByScalar(k[1], B[3][1])),
                               multiplyBsByScalar(k[2], B[3][2]))),
      dt);
  k[4] =
      multiplyDsByTime(deltaStateFunction(addBs(
                           addBs(addBs(cbs, multiplyBsByScalar(k[1], B[4][1])),
                                 multiplyBsByScalar(k[2], B[4][2])),
                           multiplyBsByScalar(k[3], B[4][3]))),
                       dt);
  k[5] = multiplyDsByTime(
      deltaStateFunction(
          addBs(addBs(addBs(addBs(cbs, multiplyBsByScalar(k[1], B[5][1])),
                            multiplyBsByScalar(k[2], B[5][2])),
                      multiplyBsByScalar(k[3], B[5][3])),
                multiplyBsByScalar(k[4], B[5][4]))),
      dt);
  k[6] = multiplyDsByTime(
      deltaStateFunction(
          addBs(addBs(addBs(addBs(addBs(cbs, multiplyBsByScalar(k[1], B[6][1])),
                                  multiplyBsByScalar(k[2], B[6][2])),
                            multiplyBsByScalar(k[3], B[6][3])),
                      multiplyBsByScalar(k[4], B[6][4])),
                multiplyBsByScalar(k[5], B[6][5]))),
      dt);
*/
  std::array<bs_struct, nOfBodies> ansBs =
      addBs(addBs(addBs(addBs(addBs(multiplyBsByScalar(k[1], CH[1]),
                                    multiplyBsByScalar(k[2], CH[2])),
                              multiplyBsByScalar(k[3], CH[3])),
                        multiplyBsByScalar(k[4], CH[4])),
                  multiplyBsByScalar(k[5], CH[5])),
            multiplyBsByScalar(k[6], CH[6]));
  cbs = addBs(cbs, ansBs);
  return cbs;
}

int main() {

  B[2][1] = 1.0 / 4.0;
  B[3][1] = 3.0 / 32;
  B[3][2] = 9.0 / 32;
  B[4][1] = 1932.0 / 2197;
  B[4][2] = -7200.0 / 2197;
  B[4][3] = 7296.0 / 2197;
  B[5][1] = 439.0 / 216;
  B[5][2] = -8;
  B[5][3] = 3680.0 / 513;
  B[5][4] = -845.0 / 4104;
  B[6][1] = -8.0 / 27;
  B[6][2] = 2;
  B[6][3] = -3544.0 / 2565;
  B[6][4] = 1859.0 / 4104;
  B[6][5] = -11.0 / 40;

  std::array<bs_struct, nOfBodies> currentBodiesState;
  currentBodiesState = g2b(currentBodiesState);
  if (hasDE) {

    SetTraceLogLevel(4);
    InitWindow(winSize.x, winSize.y, "");
    // SetTargetFPS(100);

    // setup camera

    camera.position = {camera_startPos.x, camera_startPos.y, camera_startPos.z};
    camera.target = {camera_startTarget.x, camera_startTarget.y,
                     camera_startTarget.z};
    camera.up = {camera_up.x, camera_up.y, camera_up.z};
    camera.fovy = camera_startFovy;
    camera.projection = CAMERA_PERSPECTIVE;

    Image img = GenImageColor(32, 32, WHITE);
    Texture texture = LoadTextureFromImage(img);
    UnloadImage(img);

    while (!WindowShouldClose()) {
      BeginDrawing();
      ClearBackground(BLACK);
      DrawFPS(0, 0);
      std::string txt =
          "dt: " + str(log10(constTimeStep * GetFrameTime() * speed));
      print(txt);

      BeginMode3D(camera);
      rlPushMatrix();
      rlRotatef(90, 1, 0, 0);
      drawGrid(1000, 10, WHITE);
      rlPopMatrix();
      for (int i = 0; i < nOfBodies; i++) {
        print("i=" + str(i));
        print("cbs:");
        currentBodiesState[i].Print();
        drawCircle(vecToVecf(currentBodiesState[i].position), 5, bodyColor[i],
                   drawScale);
      }
      EndMode3D();
      EndDrawing();
      currentBodiesState = RKF45(currentBodiesState);
    }

    UnloadTexture(texture);
    CloseWindow();

    std::cout << "------------WINDOW CLOSED ----------------------\n";
  }
  for (int i = 0; i < 1000; i++) {
    currentBodiesState = RKF45(currentBodiesState);
  }
  currentBodiesState[0].Print();
  std::cout << "sizeof vector: " << sizeof(vector(1, 2, 3)) << "\n";
  std::cout << "sizeof bs_struct: " << sizeof(bs_struct({1, 2}, {3, 4}, 5))
            << "\n";

  return 0;
}
