#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H
#include "vector.hpp"
#include <cmath>
#include <raylib.h>
struct {

  const double pi = 3.1415926535897932384626433832795028841971693993751;
  const double G = 4 * pow(pi, 2); // AU^3 MO^-1 Year^-2
  const double constTimeStep = pow(10, -4);
  const int circleSegments = 36;

  const Vector3 camera_startPos = {0, 0, 250};
  const Vector3 camera_startTarget = {0, 0, 0};
  const Vector3 camera_up = {0, 1, 0};
  const float camera_startFovy = 60;

  const Vector3 winSize = {1600, 900};
  const Vector3 winCenter = {winSize.x / 2, winSize.y / 2};

  const int nOfBodies = 2;
  const bool hasDE = true;

  bool trajectories = false;
  double speed = 100;

  double drawScale = 50;
  Color bodyColor[3] = {RED, BLUE, GREEN};

} global;

#endif // GLOBAL_VARIABLES_H
