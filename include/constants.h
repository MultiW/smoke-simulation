#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <Eigen/Core>

// time step
const double dt = 0.001;

// dimensions of the smoke box
const Eigen::Vector3i BOX_DIM(40, 20, 20);

// dimensions of staggered grid to compute pressure
// - NOTE: GRID_DIM - 1 (along all dimensions) must have the same 
// - propertions as BOX_DIM
const Eigen::Vector3i GRID_DIM(41, 21, 21);

// smoke particle count
const int PARTICLE_COUNT = 50;
const Eigen::AlignedBox3d SMOKE_BOUNDS(Eigen::Vector3d(5, 5, 5), Eigen::Vector3d(8, 8, 8));

// buoyancy force variables
const double FLUID_DENSITY = 2.0;
const double AIR_DENSITY = 1.0;

const double FLUID_TEMP = 1;
const double AMBIENT_TEMP = 0;

const double ALPHA = 0.05; // how much particles sink
const double BETA = 0.08; // how much particles float

//const double ALPHA = 0.08; // how much particles sink
//const double BETA = 0.97; // how much particles float

#endif // !CONSTANTS_H
