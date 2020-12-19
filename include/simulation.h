#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"
#include "grid_util.h"
#include "external_forces.h"
#include "advection.h"

#include <igl/grid.h>

#include <Eigen/Geometry>

#include <stdio.h>


// TODO: give more appropriate name for file

// Predefined colors
const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
const Eigen::RowVector3d yellow(1.0, 0.9, 0.2);
const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
const Eigen::RowVector3d green(0.2, 0.6, 0.3);
const Eigen::RowVector3d black(0.0, 0.0, 0.0);
const Eigen::RowVector3d white(1.0, 1.0, 1.0);

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
	
	// TODO: simulate pressure
}

inline void createSmokeBox(Eigen::MatrixXd& boxV, Eigen::MatrixXi& boxF, Eigen::MatrixXd& q)
{
	// Create box
	igl::read_triangle_mesh("../data/box.obj", boxV, boxF);

	Eigen::AlignedBox3d boxBounds;
	boxBounds.extend(Eigen::Vector3d(0, 0, 0));
	boxBounds.extend(Eigen::Vector3d(100, 100, 100));
	transformVertices(boxV, boxBounds);

	// Create smoke particles inside box
	// TODO: perhaps randomly create N particles within a boundary
	igl::grid(Eigen::Vector3d(20, 20, 20), q);

	Eigen::AlignedBox3d smokeBounds;
	smokeBounds.extend(Eigen::Vector3d(5, 50, 5));
	smokeBounds.extend(Eigen::Vector3d(95, 95, 95));
	transformVertices(q, smokeBounds);
}

inline void simulation_setup(int argc, char** argv, Eigen::MatrixXd& q, Eigen::MatrixXd& qdot)
{
	// Add box
	Eigen::MatrixXd boxV;
	Eigen::MatrixXi boxF;
	createSmokeBox(boxV, boxF, q);
	boxId = Visualize::addObjectToScene(boxV, boxF, orange);
	Visualize::setInvisible(boxId, true);

	// Add smoke
	smokeId = Visualize::addPointsToScene(q, white);
}

inline void draw(Eigen::Ref<const Eigen::MatrixXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, double t)
{
	Visualize::updatePoints(smokeId, q, white);
}

#endif
