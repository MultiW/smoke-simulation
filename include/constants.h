#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <Eigen/Core>

// time step
const double dt = 0.001;

// dimensions of the smoke box
const double boxX = 20;
const double boxY = 20;
const double boxZ = 20;
const Eigen::AlignedBox3d SMOKE_BOX(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(boxX, boxY, boxZ));

// dimensions of staggered grid to compute pressure
// - NOTE: GRID_DIM - 1 (along all dimensions) must have the same 
// - propertions as SMOKE_BOX
const Eigen::Vector3i GRID_DIM(21, 21, 21);

// smoke particle count
const int PARTICLE_COUNT = 50;
const Eigen::AlignedBox3d SMOKE_PARTICLE_BOUNDS(Eigen::Vector3d(2, 2, 2), Eigen::Vector3d(6, 6, 6));

// buoyancy force variables
const double FLUID_DENSITY = 2.0;
const double AIR_DENSITY = 1.0;

const double FLUID_TEMP = 1;
const double AMBIENT_TEMP = 0;

const double ALPHA = 0.05; // how much particles sink
const double BETA = 0.08; // how much particles float

//const double ALPHA = 0.08; // how much particles sink
//const double BETA = 0.97; // how much particles float

const bool ball = true;
const Eigen::Vector3d initialBallPosition(boxX/2, boxY/2, boxZ/2);
const Eigen::RowVector3d ballVelocity(100, 0, 0);
const double ballRadius = 3;

#endif // !CONSTANTS_H

