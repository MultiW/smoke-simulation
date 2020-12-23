#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <Eigen/Core>

// time step
const double dt = 0.001;

// dimensions of the smoke box
const Eigen::Vector3i BOX_DIM(200, 100, 100);

// dimensions of staggered grid to compute pressure
// - NOTE: GRID_DIM - 1 (along all dimensions) must have the same 
// - propertions as BOX_DIM
const Eigen::Vector3i GRID_DIM(21, 11, 11);

// smoke particle count
const Eigen::Vector3d SMOKE_DIM(20, 20, 20);

// buoyancy force variables
const double FLUID_DENSITY = 2.0;
const double AIR_DENSITY = 1.0;

const double FLUID_TEMP = 1;
const double AMBIENT_TEMP = 0;

const double ALPHA = 0.08;
const double BETA = 0.97;

#endif // !CONSTANTS_H
