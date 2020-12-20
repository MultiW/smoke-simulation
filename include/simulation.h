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
const Eigen::Vector3i BOX_SIZE(200, 100, 100);
// dimensions of staggered grid to compute pressure (entire vector must be factor of BOX_SIZE)
const Eigen::Vector3i GRID_SIZE(20, 10, 10);
// smoke particle count
Eigen::Vector3d SMOKE_SIZE(20, 20, 20);

// Predefined colors
const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
const Eigen::RowVector3d yellow(1.0, 0.9, 0.2);
const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
const Eigen::RowVector3d green(0.2, 0.6, 0.3);
const Eigen::RowVector3d black(0.0, 0.0, 0.0);
const Eigen::RowVector3d white(1.0, 1.0, 1.0);
const Eigen::RowVector3d red(0.8, 0.2, 0.2);

//Staggered Grid
StaggeredGrid staggeredGrid;

// Viewer data ids
int smokeId;
int boxId;

// Update location and velocity of smoke particles
inline void simulate(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt, double t)
{
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
	igl::grid(SMOKE_SIZE, q);

	Eigen::AlignedBox3d smokeBounds;
	smokeBounds.extend(Eigen::Vector3d(5, 50, 5));
	smokeBounds.extend(Eigen::Vector3d(BOX_SIZE(0) - 5, BOX_SIZE(1) - 5, BOX_SIZE(2) - 5));
	transformVertices(q, smokeBounds);
}

inline void initializeVelocity(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot)
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

inline void simulation_setup(int argc, char** argv, Eigen::MatrixXd& q, Eigen::MatrixXd& qdot)
{
	// Define boundaries of box
	Eigen::AlignedBox3d boundary;
	boundary.extend(Eigen::Vector3d(0, 0, 0));
	boundary.extend(BOX_SIZE.cast<double>());

	// Add box
	Eigen::MatrixXd boxV;
	Eigen::MatrixXi boxF;
	createSmokeBox(boxV, boxF, q, boundary);
	boxId = Visualize::addObjectToScene(boxV, boxF, orange);
	Visualize::setInvisible(boxId, true);

	// Add smoke
	smokeId = Visualize::addPointsToScene(q, white);

	// Initialize velocity
	initializeVelocity(q, qdot);

	// TODO: TEST. DELETE. START
	Eigen::MatrixXd model;
	igl::grid(GRID_SIZE, model);
	transformVertices(model, boundary);
	Visualize::addPointsToScene(model, blue);

	double cellHalfLen = boundary.sizes()(0) / (GRID_SIZE(0) - 1) / 2.0;

	Eigen::MatrixXd u, v, w, p;
	u = model;
	transformVertices(u, boundary);
	addToCol(u, 0, cellHalfLen);
	Visualize::addPointsToScene(u, yellow);

	v = model;
	transformVertices(v, boundary);
	addToCol(v, 1, cellHalfLen);
	Visualize::addPointsToScene(v, orange);

	w = model;
	transformVertices(w, boundary);
	addToCol(w, 2, cellHalfLen);
	Visualize::addPointsToScene(w, green);

	p = model;
	transformVertices(p, boundary);
	addToCol(p, 0, cellHalfLen);
	addToCol(p, 1, cellHalfLen);
	addToCol(p, 2, cellHalfLen);
	Visualize::addPointsToScene(p, red);
	// TODO: TEST. DELETE. END
}

inline void draw(Eigen::Ref<const Eigen::MatrixXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, double t)
{
	Visualize::updatePoints(smokeId, q, white);
}

#endif
