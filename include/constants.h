#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <Eigen/Core>

// time step
const double dt = 0.01;

// dimensions of the smoke box
const double boxX = 40;
const double boxY = 40;
const double boxZ = 40;
const Eigen::AlignedBox3d SMOKE_BOX(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(boxX, boxY, boxZ));

// dimensions of staggered grid to compute pressure
// - NOTE: GRID_DIM - 1 (along all dimensions) must have the same 
// - propertions as SMOKE_BOX
const Eigen::Vector3i GRID_DIM(41, 41, 41);

// smoke particle count and initial location
const int PARTICLE_COUNT = 100;
const Eigen::AlignedBox3d SMOKE_BOUNDS(Eigen::Vector3d(boxX / 2 - 3, boxY / 2 - 3, boxZ / 2 - 3), Eigen::Vector3d(boxX / 2 + 3, boxY / 2 + 3, boxZ / 2 + 3));

// buoyancy force variables
const double FLUID_DENSITY = 2.0;
const double AIR_DENSITY = 1.0;

const double FLUID_TEMP = 1;
const double AMBIENT_TEMP = 0;

const double ALPHA = 0.05; // how much particles sink
const double BETA = 2; // how much particles float

// external objects
const bool ball = false;
const Eigen::Vector3d initialBallPosition(boxX / 2, boxY / 2, boxZ / 2);
const Eigen::RowVector3d ballVelocity(20, 0, 0);
const double ballRadius = 3;

const bool bunny = false;
const Eigen::Vector3d initialBunnyPosition(boxX / 2, boxY / 2, boxZ / 2);
const Eigen::RowVector3d bunnyVelocity(20, 0, 0);
const double bunnyHalfLength = 3; // of side in the x direction

// particles
// particle center must be at (0, 0, 0)
const Eigen::AlignedBox3d PARTICLE_SIZE(Eigen::Vector3d(-0.1, -0.1, -0.1), Eigen::Vector3d(0.1, 0.1, 0.1));
const bool useParticles = true;

#endif // !CONSTANTS_H

