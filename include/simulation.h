#ifndef SIMULATION_H
#define SIMULATION_H

#include <visualization.h>
#include <igl/grid.h>

// TODO: give more appropriate name for file

// Predefined colors
const Eigen::RowVector3d orange(1.0,0.7,0.2);
const Eigen::RowVector3d yellow(1.0,0.9,0.2);
const Eigen::RowVector3d blue(0.2,0.3,0.8);
const Eigen::RowVector3d green(0.2,0.6,0.3);
const Eigen::RowVector3d black(0.0,0.0,0.0);

// Viewer data ids
int smokeBoxId;

inline void simulate(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt, double t) {

}

inline void setup(int argc, char** argv, Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) {
	// Add box
	Eigen::MatrixXd boxV;
	Eigen::MatrixXd boxF;
	igl::readOBJ("../data/box.obj", boxV, boxF);
	smokeBoxId = Visualize::addObjectToScene(boxV, boxF, orange);
	// TODO: transform size of box

	// Add smoke
	igl::grid(Eigen::Vector3d(100, 100, 100), q);
	// TODO: transform to fit inside the box
	Visualize::addPointsToScene(smokeBoxId, q, yellow);
}

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t) {

}

#endif
