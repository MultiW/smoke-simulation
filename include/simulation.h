#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"
#include "grid_util.h"
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
const Eigen::Vector3d boxSize(100, 100, 100); // NOTE: must be square
const int GRIDDIM = 20;
// Predefined colors
const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
const Eigen::RowVector3d yellow(1.0, 0.9, 0.2);
const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
const Eigen::RowVector3d green(0.2, 0.6, 0.3);
const Eigen::RowVector3d black(0.0, 0.0, 0.0);
const Eigen::RowVector3d white(1.0, 1.0, 1.0);
const Eigen::RowVector3d red(0.8, 0.2, 0.2);

//Staggered Grid
StaggeredGrid staggeredGrid(GRIDDIM);

// Viewer data ids
int smokeId;
int boxId;

// Update location and velocity of smoke particles
inline void simulate(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt, double t)
{
	// TODO: turned off simulation for now
	//// Update position q
	//advection(q, qdot, dt);

	//// Update velocity qdot
	//external_forces(q, qdot, dt);
	//
	// simulate pressure
	staggeredGrid.initVelocities(q, qdot);
	staggeredGrid.updateGrids();

}

inline void createSmokeBox(Eigen::MatrixXd& boxV, Eigen::MatrixXi& boxF, Eigen::MatrixXd& q, Eigen::AlignedBox3d& boundary)
{
	// Create box
	igl::read_triangle_mesh("../data/box.obj", boxV, boxF);

	transformVertices(boxV, boundary);

	// Create smoke particles inside box
	// TODO: perhaps randomly create N particles within a boundary
	igl::grid(Eigen::Vector3d(20, 20, 20), q);

	Eigen::AlignedBox3d smokeBounds;
	smokeBounds.extend(Eigen::Vector3d(5, 50, 5));
	smokeBounds.extend(Eigen::Vector3d(95, 95, 95));
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

// TODO: CREATE NEW UTIL FUNC
inline void addToCol(Eigen::MatrixXd& matrix, int columnIndex, double value)
{
	Eigen::VectorXd ones;
	ones.resize(matrix.rows(), 1);
	ones.setOnes();
	matrix.col(columnIndex) += ones * value;
}

inline void simulation_setup(int argc, char** argv, Eigen::MatrixXd& q, Eigen::MatrixXd& qdot)
{
	// Define boundaries of box
	Eigen::AlignedBox3d boundary;
	boundary.extend(Eigen::Vector3d(0, 0, 0));
	boundary.extend(boxSize);

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
	Eigen::MatrixXd g;
	igl::grid(Eigen::Vector3d(20, 20, 20), g);
	Eigen::AlignedBox3d gBox;
	createAlignedBox(g, gBox);
	transformVertices(g, boundary);
	Visualize::addPointsToScene(g, blue);


	double cellHalfLen = boundary.sizes()(0) / 19.0 / 2.0;

	Eigen::MatrixXd u, v, w, p;

	igl::grid(Eigen::Vector3d(GRIDDIM - 1, GRIDDIM, GRIDDIM), u);
	transformVertices(u, gBox, boundary, false);
	addToCol(u, 0, cellHalfLen);
	Visualize::addPointsToScene(u, yellow);

	igl::grid(Eigen::Vector3d(GRIDDIM, GRIDDIM - 1, GRIDDIM), v);
	transformVertices(v, gBox, boundary, false);
	addToCol(v, 1, cellHalfLen);
	Visualize::addPointsToScene(v, orange);

	igl::grid(Eigen::Vector3d(GRIDDIM, GRIDDIM, GRIDDIM - 1), w);
	transformVertices(w, gBox, boundary, false);
	addToCol(w, 2, cellHalfLen);
	Visualize::addPointsToScene(w, green);

	igl::grid(Eigen::Vector3d(GRIDDIM - 1, GRIDDIM - 1, GRIDDIM - 1), p);
	transformVertices(p, gBox, boundary, false);
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
