#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"

#include "grid_util.h"
#include "util.h"
#include "staggered_grid.h"
#include "constants.h"
#include "util.h"

#include <igl/grid.h>
#include <igl/bounding_box.h>

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
int bunnyId;

// == Simulation State ==
Eigen::MatrixXd q;

Eigen::RowVector3d currBallCenter;
Eigen::MatrixXd ballV;
Eigen::MatrixXi ballF;

Eigen::RowVector3d currBunnyCenter;
Eigen::MatrixXd bunnyV;
Eigen::MatrixXi bunnyF;
// ====

// Update location and velocity of smoke particles
inline void simulate()
{
	staggeredGrid.updateExternalObjects(currBallCenter, &bunnyV, &bunnyF);

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

inline void createBall(Eigen::MatrixXd& ballV, Eigen::MatrixXi& ballF)
{
   	igl::read_triangle_mesh("../data/sphere.obj", ballV, ballF);

   	// Display sphere smaller than actual size to account for the particle's large size
   	Eigen::AlignedBox3d sphereBoundaries;
   	sphereBoundaries.extend(initialBallPosition - Eigen::Vector3d::Constant(ballRadius - 0.5));
   	sphereBoundaries.extend(initialBallPosition + Eigen::Vector3d::Constant(ballRadius - 0.5));
   	transformVertices(ballV, sphereBoundaries);
}

inline void createBunny(Eigen::MatrixXd& bunnyV, Eigen::MatrixXi bunnyF)
{
	igl::read_triangle_mesh("../data/bunny.off", bunnyV, bunnyF);

	// Find dimensions of bunny
	Eigen::MatrixXd boundingV;
	Eigen::MatrixXi boundingF;
	igl::bounding_box(bunnyV, boundingV, boundingF);
	Eigen::AlignedBox3d defaultBox;
	createAlignedBox(boundingV, defaultBox);

	// Scale bunny to set size
	double scale = bunnyHalfLengths / (defaultBox.sizes()(0) / 2.0);
	double xHalfLen = bunnyHalfLengths;
	double yHalfLen = defaultBox.sizes()(1) * scale / 2;
	double zHalfLen = defaultBox.sizes()(2) * scale / 2;
	//Eigen::Vector3d minCorner = initialBunnyPosition - ;
	Eigen::Vector3d maxCorner;
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
		currBallCenter = initialBallPosition.transpose();
		createBall(ballV, ballF);
		ballId = Visualize::addObjectToScene(ballV, ballF, orange);
	}
	if (bunny)
	{
		currBunnyCenter = initialBunnyPosition.transpose();
		createBunny(bunnyV, bunnyF);
		bunnyId = Visualize::addObjectToScene(bunnyV, bunnyF, orange);
	}

	Visualize::viewer().selected_data_index = Visualize::viewer().mesh_index(boxId);

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
