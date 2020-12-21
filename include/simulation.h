#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"

#include "grid_util.h"
#include "util.h"
#include "external_forces.h"
#include "advection.h"
#include "staggered_grid.h"

#include <igl/grid.h>

#include <Eigen/Geometry>

#include <stdio.h>
#include <time.h>
#include <iostream>


// TODO: give more appropriate name for file

// Tuning parameters
// dimensions of the smoke box
const Eigen::Vector3i BOX_DIM(200, 100, 100);
// dimensions of staggered grid to compute pressure
// - NOTE: GRID_DIM - 1 (along all dimensions) must have the same 
// - propertions as BOX_DIM
const Eigen::Vector3i GRID_DIM(21, 11, 11);
// smoke particle count
Eigen::Vector3d SMOKE_DIM(20, 20, 20);

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
inline void simulate(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt, double t)
{
	// TODO: add boundary checks on points

	// Update position q
	advection(q, qdot, dt);

	// Update velocity qdot
	external_forces(q, qdot, dt);

	// simulate pressure
	staggeredGrid.computeVelocity(q, qdot);
}

inline void createSmokeBox(Eigen::MatrixXd& boxV, Eigen::MatrixXi& boxF, Eigen::MatrixXd& q, Eigen::AlignedBox3d& boundary)
{
	// Create box
	igl::read_triangle_mesh("../data/box.obj", boxV, boxF);

	transformVertices(boxV, boundary);

	// Create smoke particles inside box
	// TODO: perhaps randomly create N particles within a boundary
	igl::grid(SMOKE_DIM, q);

	Eigen::AlignedBox3d smokeBounds;
	smokeBounds.extend(Eigen::Vector3d(5, 50, 5));
	smokeBounds.extend(Eigen::Vector3d(BOX_DIM(0) - 5, BOX_DIM(1) - 5, BOX_DIM(2) - 5));
	transformVertices(q, smokeBounds);
}

inline void initParticleVelocity(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot)
{
	qdot.resize(q.rows(), q.cols());
	std::srand((unsigned) std::time(NULL));
	for (int j = 0; j < qdot.cols(); j++)
	{
		for (int i = 0; i < qdot.rows(); i++)
		{
			// iterate in column-major order, the default order for Eigen matrix
			qdot(i, j) = (double) std::rand() / RAND_MAX * 100.0;
		}
	}
}

// Must be called first
inline void simulation_setup(int argc, char** argv, Eigen::MatrixXd& q, Eigen::MatrixXd& qdot)
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

	initParticleVelocity(q, qdot);

	staggeredGrid = StaggeredGrid(smokeBox, GRID_DIM);

	// TODO: DELETE. Testing if initialization of staggered grid points is correct
	Eigen::MatrixXd u, v, w, p;
	staggeredGrid.getGridPoints(u, v, w, p);
	Visualize::addPointsToScene(u, yellow);
	Visualize::addPointsToScene(v, orange);
	Visualize::addPointsToScene(w, green);
	Visualize::addPointsToScene(p, red);
}

inline void draw(Eigen::Ref<const Eigen::MatrixXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, double t)
{
	Visualize::updatePoints(smokeId, q, white);
}

#endif
