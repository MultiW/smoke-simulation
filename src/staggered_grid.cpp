#include "staggered_grid.h"





StaggeredGrid::StaggeredGrid(int dim): dim(dim) {
	this->uGrid.resize(dim - 1, dim, dim);
	this->vGrid.resize(dim, dim - 1, dim);
	this->wGrid.resize(dim, dim, dim - 1);
	this->pGrid.resize(dim - 1, dim - 1, dim - 1);
}

void StaggeredGrid::initVelocities(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot) {
	//TODO interpolate qdot into grids
}

void StaggeredGrid::updateGrids() {

}