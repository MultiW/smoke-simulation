#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"

#include "grid_util.h"
#include "util.h"
#include "external_forces.h"
#include "advection.h"
#include "staggered_grid.h"
#include "constants.h"

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
inline void simulate(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double t)
{
	// 1. update velocities
	staggeredGrid.advectVelocities();
	staggeredGrid.applyExternalForces();
	staggeredGrid.applyPressureProjections();

	// 2. advect temperature and density
	staggeredGrid.updateTemperatureAndDensity();

	// 3. advect particles
	// TODO: move to staggered_grid
	//  - 1. find v^t: current velocity of particle (using Grid.interpolate)
	//  - 2. find v^t+1: find particle's next position (+ dt * v^t), then find velocity at that position (using Grid.interpolate)
	//  - 3. average v^t and v^t+1
	//  - 4. compute next position of particle using: + dt * v_avg
	//  - Whenever a new position is computed: clip position of particle to the grid
	advection(q, qdot, dt);
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
	smokeBounds.extend(Eigen::Vector3d(5, BOX_DIM(1) - 20, 5));
	smokeBounds.extend(Eigen::Vector3d(BOX_DIM(0) - 5, BOX_DIM(1) - 5, BOX_DIM(2) - 5));
	transformVertices(q, smokeBounds);
}

inline void initParticleVelocity(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot)
{
	// TODO: how should velocity be initialized
	qdot.resize(q.rows(), q.cols());
	qdot.setZero();
	qdot.col(1).setConstant(-0.01);
	//std::srand((unsigned) std::time(NULL));
	//for (int j = 0; j < qdot.cols(); j++)
	//{
	//	for (int i = 0; i < qdot.rows(); i++)
	//	{
	//		// iterate in column-major order, the default order for Eigen matrix
	//		qdot(i, j) = (double) std::rand() / RAND_MAX * 100.0;
	//	}
	//}
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

	staggeredGrid = StaggeredGrid(smokeBox, GRID_DIM);

	//// TODO: DELETE. Testing if initialization of staggered grid points is correct
	//Eigen::MatrixXd u, v, w, p;
	//staggeredGrid.getGridPoints(u, v, w, p);
	//Visualize::addPointsToScene(u, yellow);
	//Visualize::addPointsToScene(v, orange);
	//Visualize::addPointsToScene(w, green);
	//Visualize::addPointsToScene(p, red);
}

inline void draw(Eigen::Ref<const Eigen::MatrixXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, double t)
{
	Visualize::updatePoints(smokeId, q, white);
}

#endif
