#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"

#include "grid_util.h"
#include "util.h"
#include "external_forces.h"
#include "staggered_grid.h"
#include "constants.h"
#include "util.h"

#include <igl/grid.h>

#include <Eigen/Geometry>

#include <stdio.h>
#include <time.h>
#include <iostream>


// Predefined colors
const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
const Eigen::RowVector3d yellow(1.0, 0.9, 0.2);
const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
const Eigen::RowVector3d green(0.2, 0.6, 0.3);
const Eigen::RowVector3d black(0.0, 0.0, 0.0);
const Eigen::RowVector3d white(1.0, 1.0, 1.0);
const Eigen::RowVector3d red(0.8, 0.2, 0.2);

// Helper variables
StaggeredGrid staggeredGrid;

// Viewer data ids
int smokeId;
int boxId;

// Update location and velocity of smoke particles
inline void simulate(Eigen::MatrixXd& q, double t)
{
	// 1. update velocities
	staggeredGrid.advectVelocities();
	staggeredGrid.applyExternalForces();
	staggeredGrid.applyPressureProjections();

	// 2. advect temperature and density
	staggeredGrid.updateTemperatureAndDensity();

	// 3. advect particles
	staggeredGrid.advectPosition(q);
}

inline void createSmokeBox(Eigen::MatrixXd& boxV, Eigen::MatrixXi& boxF, Eigen::MatrixXd& q, Eigen::AlignedBox3d& boundary)
{
	// Create box
	igl::read_triangle_mesh("../data/box.obj", boxV, boxF);

	transformVertices(boxV, boundary);

	// Create smoke particles inside box
	q.resize(PARTICLE_COUNT, 3);
	Eigen::Vector3d bottomLeftFloor = SMOKE_BOUNDS.corner(Eigen::AlignedBox3d::BottomLeftFloor);
	Eigen::Vector3d topRightCeil = SMOKE_BOUNDS.corner(Eigen::AlignedBox3d::TopRightCeil);
	for (int i = 0; i < PARTICLE_COUNT; i++)
	{
		q(i, 0) = getRand(bottomLeftFloor(0), topRightCeil(0));
		q(i, 1) = getRand(bottomLeftFloor(1), topRightCeil(1));
		q(i, 2) = getRand(bottomLeftFloor(2), topRightCeil(2));
	}
}

// Must be called first
inline void simulation_setup(int argc, char** argv, Eigen::MatrixXd& q)
{
	// Define boundaries of box
	Eigen::AlignedBox3d smokeBox;
	smokeBox.extend(Eigen::Vector3d(0, 0, 0));
	smokeBox.extend(BOX_DIM.cast<double>());

	// Add box
	Eigen::MatrixXd boxV;
	Eigen::MatrixXi boxF;
	createSmokeBox(boxV, boxF, q, smokeBox);
	boxId = Visualize::addObjectToScene(boxV, boxF, orange);
	Visualize::setInvisible(boxId, true);

	smokeId = Visualize::addPointsToScene(q, white);

	staggeredGrid = StaggeredGrid(smokeBox, GRID_DIM);

	////// TODO: DELETE. Testing if initialization of staggered grid points is correct
	//Eigen::MatrixXd u, v, w, p;
	//staggeredGrid.getGridPoints(u, v, w, p);
	//Visualize::addPointsToScene(u, yellow);
	//Visualize::addPointsToScene(v, orange);
	//Visualize::addPointsToScene(w, green);
	//Visualize::addPointsToScene(p, red);
}

inline void draw(Eigen::Ref<const Eigen::MatrixXd> q, double t)
{
	Visualize::updatePoints(smokeId, q, white);
}

#endif
