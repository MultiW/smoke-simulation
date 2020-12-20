#include "staggered_grid.h"

StaggeredGrid::StaggeredGrid(size_t dim) :
	dim(dim),
	uGrid(Grid(dim - 1, dim, dim)),
	vGrid(Grid(dim, dim - 1, dim)),
	wGrid(Grid(dim, dim, dim - 1)),
	pGrid(Grid(dim - 1, dim - 1, dim - 1)) {

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

