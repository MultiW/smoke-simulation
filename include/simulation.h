#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"
#include "grid_util.h"

#include <igl/grid.h>

#include <Eigen/Geometry>


// TODO: give more appropriate name for file

// Predefined colors
const Eigen::RowVector3d orange(1.0,0.7,0.2);
const Eigen::RowVector3d yellow(1.0,0.9,0.2);
const Eigen::RowVector3d blue(0.2,0.3,0.8);
const Eigen::RowVector3d green(0.2,0.6,0.3);
const Eigen::RowVector3d black(0.0,0.0,0.0);

// Viewer data ids
int smokeBoxId;

inline void simulate(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt, double t) 
{

}

inline void createBox(Eigen::MatrixXd& boxV, Eigen::MatrixXi& boxF)
{
//	boxV.resize(8, 3);
//	boxV <<
//		0, 0, 0,
//		0, 0, 1,
//		0, 1, 0,
//		0, 1, 1,
//		1, 0, 0,
//		1, 0, 1,
//		1, 1, 0,
//		1, 1, 1;
//	boxF.resize(12, 3);
//	boxF <<
//		0, 1, 3,
//		3, 1, 2,
////		0, 1, 2,
////		2, 1, 3;

	smokeBoxId = Visualize::addObjectToScene(boxV, boxF, orange);

	Eigen::AlignedBox3d bounds;
	bounds.extend(Eigen::Vector3d(0, 0, 0));
	bounds.extend(Eigen::Vector3d(100, 100, 100));
}

inline void simulation_setup(int argc, char** argv, Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) 
{
	// Add box
	Eigen::MatrixXd boxV;
	Eigen::MatrixXi boxF;
	createBox(boxV, boxF);

	// Add smoke
	igl::grid(Eigen::Vector3d(100, 100, 100), q);
	// TODO: transform to fit inside the box
	Visualize::addPointsToScene(smokeBoxId, q, yellow);
}

inline void draw(Eigen::Ref<const Eigen::MatrixXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, double t) {
	// Necessary? Think it's only useful for skinning
}

#endif
