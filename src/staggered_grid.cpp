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
	uGrid(Grid(dim(0) - 1.0, dim(1), dim(2))),
	vGrid(Grid(dim(0), dim(1) - 1.0, dim(2))),
	wGrid(Grid(dim(0), dim(1), dim(2) - 1.0)),
	pGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0))
{
	// World-space locations of grid-points
	Eigen::MatrixXd u, v, w, p;
	this->createGridPoints(u, v, w, p);

	// Update world points of grids
	this->updateWorldPoints(u, this->uGrid);
	this->updateWorldPoints(v, this->vGrid);
	this->updateWorldPoints(w, this->wGrid);
	this->updateWorldPoints(p, this->pGrid);
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
	addToCol(u, 0, cellHalfLen);

	v = model;
	transformVertices(v, this->box);
	addToCol(v, 1, cellHalfLen);

	w = model;
	transformVertices(w, this->box);
	addToCol(w, 2, cellHalfLen);

	p = model;
	transformVertices(p, this->box);
	addToCol(p, 0, cellHalfLen);
	addToCol(p, 1, cellHalfLen);
	addToCol(p, 2, cellHalfLen);
}

int getPointIdx(std::vector<double> &sortedPoints, double interval, double currPoint, double epsilon)
{
	double min = sortedPoints.front();
	double max = sortedPoints.back();

	// Check that point is in bounds
	assert(currPoint >= min - epsilon && currPoint <= max + epsilon);

	int borderIdx = std::round((currPoint - min) / interval);

	// Check that point is within "epsilon" from nearest grid point
	assert(std::abs(currPoint - sortedPoints[borderIdx]) < epsilon);

	return borderIdx;
}

void StaggeredGrid::updateWorldPoints(const Eigen::MatrixXd& points, Grid& grid)
{
	// Grid's sorted x, y, z values
	// Note: Eigen Matrix is stored in column-major order
	std::set<double> xSet(points.col(0).data(),	points.col(0).data() + points.rows());
	std::set<double> ySet(points.col(1).data(),	points.col(1).data() + points.rows());
	std::set<double> zSet(points.col(2).data(),	points.col(2).data() + points.rows());
	std::vector<double> xPoints(xSet.begin(), xSet.end());
	std::vector<double> yPoints(ySet.begin(), ySet.end());
	std::vector<double> zPoints(zSet.begin(), zSet.end());
	std::sort(xPoints.begin(), xPoints.end());
	std::sort(yPoints.begin(), yPoints.end());
	std::sort(zPoints.begin(), zPoints.end());

	double cellSize = this->getCellSize();
	double epsilon = cellSize / 10.0;

	int xIdx, yIdx, zIdx;
	Eigen::RowVector3d currPoint;
	for (int i = 0; i < points.rows(); i++)
	{
		currPoint = points.row(i);

		xIdx = getPointIdx(xPoints, cellSize, currPoint(0), epsilon);
		yIdx = getPointIdx(yPoints, cellSize, currPoint(1), epsilon);
		zIdx = getPointIdx(zPoints, cellSize, currPoint(2), epsilon);

		if (xIdx < grid.getSize(0) && yIdx < grid.getSize(1) && zIdx < grid.getSize(2))
		{
			grid(xIdx, yIdx, zIdx).gridPoint = Eigen::Vector3i(xIdx, yIdx, zIdx);
			grid(xIdx, yIdx, zIdx).worldPoint = Eigen::Vector3d(currPoint(0), currPoint(1), currPoint(2));
		}
		else
		{
			// TODO: REMOVE. For debugging
			printf("Not adding point: %d %d %d\n", xIdx, yIdx, zIdx);
		}
	}
}

void convertGridToPoints(Grid grid, Eigen::MatrixXd& points)
{
	points.resize(0, 3);
	for (int x = 0; x < grid.getSize(0); x++)
	{
		for (int y = 0; y < grid.getSize(1); y++)
		{
			for (int z = 0; z < grid.getSize(2); z++)
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

