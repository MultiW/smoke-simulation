#include "staggered_grid.h"

StaggeredGrid::StaggeredGrid() {}

StaggeredGrid::StaggeredGrid(Eigen::Vector3i dim) :
	dim(dim),
	uGrid(Grid(dim(0) - 1, dim(1), dim(2))),
	vGrid(Grid(dim(0), dim(1) - 1, dim(2))),
	wGrid(Grid(dim(0), dim(1), dim(2) - 1)),
	pGrid(Grid(dim(0) - 1, dim(1) - 1, dim(2) - 1))
{

}

void StaggeredGrid::setGridVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) {
	//TODO interpolate qdot into grids
}

void StaggeredGrid::getVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) {
	// TODO: trillinear interpolation
}


void StaggeredGrid::updateVelocityAndPressure() {

}

void StaggeredGrid::computeVelocity(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) {
	this->setGridVelocities(q, qdot);
	this->updateVelocityAndPressure();
	this->getVelocities(q, qdot);
}

