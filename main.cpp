#include "global_variables.hpp"
#include "vector.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <raylib.h>
#include <rlgl.h>
#include <string>

double B[7][6];

double CH[] = {0,         16.0 / 135, 0, 6656.0 / 12825, 28561.0 / 56430,
               -9.0 / 50, 2.0 / 55};
void print(std::string out) { std::cout << out << "\n"; }
std::string str(double num) { return std::to_string(num); }

struct bs_struct {
  mml::vector2 position;
  mml::vector2 velocity;
  double mass;

  bs_struct(mml::vector2 p, mml::vector2 v, double m)
      : position(p), velocity(v), mass(m) {}
  bs_struct() : position({0, 0}), velocity({0, 0}), mass(0){};
  void Print() {
    print("PRINT bs_struct: START--------------------");
    print("position:");
    printvec(position);
    print("velocity:");
    printvec(velocity);
    print("mass: " + str(mass));
    print("PRINT bs_struct: END----------------------");
  }
};

struct ds_struct {
  mml::vector2 velocity;
  mml::vector2 acceleration;

  ds_struct(mml::vector2 v, mml::vector2 a) : velocity(v), acceleration(a) {}
  ds_struct() : velocity({0, 0}), acceleration({0, 0}){};

  void Print() {
    print("PRINT ds_struct: START--------------------");
    print("velocity:");
    printvec(velocity);
    print("acceleration:");
    printvec(acceleration);
    print("PRINT ds_struct: END----------------------");
  }
};

const int nOfBodies = 2;
const double pi = global.pi;
const double G = global.G;
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
/*
void drawCircleOld(vector centerOffset, double radius, Color color) {
  int segments = global.circleSegments;
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
}*/
void drawCircle(mml::vector2 centerOffset, float radius, Color color,
                float scale) {
  int segments = global.circleSegments;
  double step = 2 * pi / segments;
  double angle = 0;
  mml::vector2 center = {centerOffset.x * scale, centerOffset.y * scale};
  mml::vector2 vertex2 = {sinf(angle + step) * (float)radius,
                          cosf(angle + step) * (float)radius};
  mml::vector2 vertex2old = {0, (float)radius};
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

std::array<ds_struct, nOfBodies>
deltaStateFunction(std::array<bs_struct, nOfBodies> bsLoc) {
  mml::vector2 ForceArr[nOfBodies][nOfBodies];
  mml::vector2 Force[nOfBodies];
  for (int i = 0; i < nOfBodies; i++) {
    for (int j = 0; j < nOfBodies; j++) {
      if (i == j)
        continue;
      if (i > j) {
        ForceArr[i][j] = ForceArr[j][i].opposite();
        continue;
      }
      mml::vector2 distance = bsLoc[j].position - bsLoc[i].position;
      double distanceMagnitude = distance.mag();
      mml::vector2 direction = distance.unit();
      double forceMagnitude = 0;
      if (distanceMagnitude > 0) {
        forceMagnitude =
            (G * bsLoc[i].mass * bsLoc[j].mass) / pow(distanceMagnitude, 2);
      }
      ForceArr[i][j] = direction * forceMagnitude;
    }
    for (int j = 0; j < nOfBodies; j++) {
      if (i == j)
        continue;
      Force[i] = Force[i] + ForceArr[i][j];
    }
  }
  std::array<ds_struct, nOfBodies> out;
  for (int i = 0; i < nOfBodies; i++) {
    out[i] = ds_struct(bsLoc[i].velocity, Force[i] / bsLoc[i].mass);
  }
  return out;
}

std::array<bs_struct, nOfBodies>
multiplyDsByTime(std::array<ds_struct, nOfBodies> ds, double t) {
  std::array<bs_struct, nOfBodies> out;
  for (int i = 0; i < nOfBodies; i++) {
    out[i] = bs_struct(ds[i].velocity * t, ds[i].acceleration * t, 0);
  }
  return out;
}
std::array<bs_struct, nOfBodies>
multiplyBsByScalar(std::array<bs_struct, nOfBodies> bs, double scalar) {
  std::array<bs_struct, nOfBodies> out;
  for (int i = 0; i < nOfBodies; i++) {
    out[i] = bs_struct(bs[i].position * scalar, bs[i].velocity * scalar,
                       bs[i].mass * scalar);
  }
  return out;
}
std::array<bs_struct, nOfBodies> addBs(std::array<bs_struct, nOfBodies> bs1,
                                       std::array<bs_struct, nOfBodies> bs2) {
  std::array<bs_struct, nOfBodies> out;
  for (int i = 0; i < nOfBodies; i++) {
    out[i] =
        bs_struct(bs1[i].position + bs2[i].position,
                  bs1[i].velocity + bs2[i].velocity, bs1[i].mass + bs2[i].mass);
  }
  return out;
}
std::array<bs_struct, nOfBodies> RKF45(std::array<bs_struct, nOfBodies> cbs) {
  double dt = global.constTimeStep * GetFrameTime() * global.speed;
  if (!global.hasDE) {
    dt = global.constTimeStep;
  }
  std::array<std::array<bs_struct, nOfBodies>, 7> k;

  for (int i = 1; i < 7; i++) {
    std::array<bs_struct, nOfBodies> dsfArg = cbs;
    for (int j = 1; j < i; j++) {
      dsfArg = addBs(dsfArg, multiplyBsByScalar(k[j], B[i][j]));
    }
    k[i] = multiplyDsByTime(deltaStateFunction(dsfArg), dt);
  }
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
  if (global.hasDE) {

    SetTraceLogLevel(4);
    InitWindow(global.winSize.x, global.winSize.y, "");
    // SetTargetFPS(100);
    // setup camera

    camera.position = global.camera_startPos;
    camera.target = global.camera_startTarget;
    camera.up = global.camera_up;
    camera.fovy = global.camera_startFovy;
    camera.projection = CAMERA_PERSPECTIVE;

    Image img = GenImageColor(32, 32, WHITE);
    Texture texture = LoadTextureFromImage(img);
    UnloadImage(img);

    while (!WindowShouldClose()) {
      BeginDrawing();
      ClearBackground(BLACK);
      DrawFPS(0, 0);
      std::string txt =
          "log10(dt): " +
          str(log10(global.constTimeStep * GetFrameTime() * global.speed));
      print(txt);

      BeginMode3D(camera);
      rlPushMatrix();
      rlRotatef(90, 1, 0, 0);
      drawGrid(1000, 10, WHITE);
      rlPopMatrix();
      for (int i = 0; i < nOfBodies; i++) {
        drawCircle(currentBodiesState[i].position, 5, global.bodyColor[i],
                   global.drawScale);
      }
      EndMode3D();
      EndDrawing();
      currentBodiesState = RKF45(currentBodiesState);
    }

    UnloadTexture(texture);
    CloseWindow();

    std::cout << "------------WINDOW CLOSED ----------------------\n";
  }
  for (int i = 0; i < 1000000; i++) {
    currentBodiesState = RKF45(currentBodiesState);
  }
  currentBodiesState[0].Print();
  std::cout << "sizeof vector: " << sizeof(mml::vector2(1, 2)) << "\n";
  std::cout << "sizeof bs_struct: " << sizeof(bs_struct({1, 2}, {3, 4}, 5))
            << "\n";

  return 0;
}
