#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"

#include "grid_util.h"
#include "util.h"
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
int ballId;

// Simulation state
Eigen::MatrixXd q;
Eigen::MatrixXd ballV;
Eigen::MatrixXi ballF;
Eigen::RowVector3d currBallCenter;

// Update location and velocity of smoke particles
inline void simulate()
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

inline void simulateBall()
{
	Eigen::RowVector3d ballMaxCorner = currBallCenter + Eigen::RowVector3d::Constant(ballRadius) + dt * ballVelocity;
	Eigen::RowVector3d ballMinCorner = currBallCenter - Eigen::RowVector3d::Constant(ballRadius) + dt * ballVelocity;
	if (isInBox(SMOKE_BOX, ballMinCorner) && isInBox(SMOKE_BOX, ballMaxCorner))
	{
		for (int i = 0; i < ballV.rows(); i++)
		{
			ballV.row(i) += dt * ballVelocity;
		}
		currBallCenter += dt * ballVelocity;
	}
}

inline void createSmokeBox(Eigen::MatrixXd& boxV, Eigen::MatrixXi& boxF, Eigen::MatrixXd& q, const Eigen::AlignedBox3d& boundary)
{
	// Create box
	igl::read_triangle_mesh("../data/box.obj", boxV, boxF);

	transformVertices(boxV, boundary);

	// Create smoke particles inside box
	q.resize(PARTICLE_COUNT, 3);
	Eigen::Vector3d bottomLeftFloor = SMOKE_PARTICLE_BOUNDS.corner(Eigen::AlignedBox3d::BottomLeftFloor);
	Eigen::Vector3d topRightCeil = SMOKE_PARTICLE_BOUNDS.corner(Eigen::AlignedBox3d::TopRightCeil);
	for (int i = 0; i < PARTICLE_COUNT; i++)
	{
		q(i, 0) = getRand(bottomLeftFloor(0), topRightCeil(0));
		q(i, 1) = getRand(bottomLeftFloor(1), topRightCeil(1));
		q(i, 2) = getRand(bottomLeftFloor(2), topRightCeil(2));
	}
}

// Must be called first
inline void simulation_setup(int argc, char** argv)
{
	// Add box
	Eigen::MatrixXd boxV;
	Eigen::MatrixXi boxF;
	createSmokeBox(boxV, boxF, q, SMOKE_BOX);
	boxId = Visualize::addObjectToScene(boxV, boxF, orange);
	Visualize::setInvisible(boxId, true);

	smokeId = Visualize::addPointsToScene(q, white);

	staggeredGrid = StaggeredGrid(q, SMOKE_BOX, GRID_DIM, SMOKE_BOX.sizes()(0) / (GRID_DIM(0) - 1.0));

	// Add sphere
	if (ball)
	{
		currBallCenter = initialBallPosition;
		igl::read_triangle_mesh("../data/sphere.obj", ballV, ballF);

		Eigen::AlignedBox3d sphereBoundaries;
		sphereBoundaries.extend(initialBallPosition - Eigen::Vector3d::Constant(ballRadius));
		sphereBoundaries.extend(initialBallPosition + Eigen::Vector3d::Constant(ballRadius));
		transformVertices(ballV, sphereBoundaries);

		ballId = Visualize::addObjectToScene(ballV, ballF, orange);
	}

	Visualize::viewer().core().align_camera_center(ballV);

	////// TODO: DELETE. Testing if initialization of staggered grid points is correct
	//Eigen::MatrixXd u, v, w, p;
	//staggeredGrid.getGridPoints(u, v, w, p);
	//Visualize::addPointsToScene(u, yellow);
	//Visualize::addPointsToScene(v, orange);
	//Visualize::addPointsToScene(w, green);
	//Visualize::addPointsToScene(p, red);
}

inline void draw()
{
	Visualize::updatePoints(smokeId, q, white);
	Visualize::updateObject(ballId, ballV);
}

#endif
