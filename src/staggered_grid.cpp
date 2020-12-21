#include "staggered_grid.h"
#include "util.h"
#include "grid_util.h"

#include <igl/grid.h>

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <set>


StaggeredGrid::StaggeredGrid() {}

StaggeredGrid::StaggeredGrid(const Eigen::AlignedBox3d& box, const Eigen::Vector3i& dim) :
	box(box),
	dim(dim),
	uGrid(Grid(dim(0), dim(1) - 1.0, dim(2) - 1.0)),
	vGrid(Grid(dim(0)- 1.0, dim(1), dim(2) - 1.0)),
	wGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2))),
	pGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0))
{
	// World-space locations of grid-points
	Eigen::MatrixXd u, v, w, p;
	this->createGridPoints(u, v, w, p);

	// Update world points of grids
	this->uGrid.setWorldPoints(u, this->getCellSize());
	this->vGrid.setWorldPoints(v, this->getCellSize());
	this->wGrid.setWorldPoints(w, this->getCellSize());
	this->pGrid.setWorldPoints(p, this->getCellSize());
}

double StaggeredGrid::getCellSize()
{
	return this->box.sizes()(0) / (this->dim(0) - 1.0);
}

// =====================================
// === Setting world point locations ===
// =====================================

void StaggeredGrid::createGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p)
{
	Eigen::MatrixXd model;
	igl::grid(this->dim, model);
	transformVertices(model, this->box);

	double cellHalfLen = this->getCellSize() / 2.0;

	u = model;
	transformVertices(u, this->box);
	addToCol(u, 1, cellHalfLen);
	addToCol(u, 2, cellHalfLen);

	v = model;
	transformVertices(v, this->box);
	addToCol(v, 0, cellHalfLen);
	addToCol(v, 2, cellHalfLen);

	w = model;
	transformVertices(w, this->box);
	addToCol(w, 0, cellHalfLen);
	addToCol(w, 1, cellHalfLen);

	p = model;
	transformVertices(p, this->box);
	addToCol(p, 0, cellHalfLen);
	addToCol(p, 1, cellHalfLen);
	addToCol(p, 2, cellHalfLen);
}

void convertGridToPoints(Grid grid, Eigen::MatrixXd& points)
{
	points.resize(0, 3);
	for (int x = 0; x < grid.size(0); x++)
	{
		for (int y = 0; y < grid.size(1); y++)
		{
			for (int z = 0; z < grid.size(2); z++)
			{
				addRows(points, grid(x, y, z).worldPoint.transpose());
			}
		}
	}
}

void StaggeredGrid::getGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p)
{
	convertGridToPoints(this->uGrid, u);
	convertGridToPoints(this->vGrid, v);
	convertGridToPoints(this->wGrid, w);
	convertGridToPoints(this->pGrid, p);
}


// ======================================
// === Computing pressure projections ===
// ======================================

void StaggeredGrid::setGridVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) {
	//TODO interpolate qdot into grids
}

void StaggeredGrid::getVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) 
{
	this->uGrid.interpolateToPoints(q, qdot.col(0));
	this->vGrid.interpolateToPoints(q, qdot.col(1));
	this->wGrid.interpolateToPoints(q, qdot.col(2));
}


void StaggeredGrid::updateVelocityAndPressure() {

}

void StaggeredGrid::computeVelocity(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) {
	// TODO: uncomment when all have been implemented
	//this->setGridVelocities(q, qdot);
	//this->updateVelocityAndPressure();
	//this->getVelocities(q, qdot);
}

