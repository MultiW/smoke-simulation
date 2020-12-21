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
	std::vector<double> xPoints;
	std::vector<double> yPoints;
	std::vector<double> zPoints;
	colToSortedVector(points.col(0), xPoints);
	colToSortedVector(points.col(1), yPoints);
	colToSortedVector(points.col(2), zPoints);

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

		if (xIdx < grid.size(0) && yIdx < grid.size(1) && zIdx < grid.size(2))
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

// binBorders - sorted array of boundaries of all bins
int getBinIdx(std::vector<double> binBorders, double binSize, double currLocation)
{
	int min = binBorders.front();
	int max = binBorders.back();
	assert(currLocation >= min && currLocation <= max);

	for (int i = 1; i < binBorders.size(); i++)
	{
		if (currLocation < binBorders[i])
		{
			return i - 1;
		}
	}

	printf("Error in getBinIdx(): could not find bin index.");
	throw;
}

void StaggeredGrid::interpolateGrid(Eigen::MatrixXd& q, Eigen::VectorXd& qdotCol, Grid& grid)
{
	double cellLen = this->getCellSize();

	// TODO: fix this. vector is wrong
	std::vector<double> xPoints;
	std::vector<double> yPoints;
	std::vector<double> zPoints;
	colToSortedVector(q.col(0), xPoints);
	colToSortedVector(q.col(1), yPoints);
	colToSortedVector(q.col(2), zPoints);

	Eigen::RowVector3d point;
	int xIdx, yIdx, zIdx; // cell in which point is located
	for (int i = 0; i < q.rows(); i++)
	{
		point = q.row(i);

		// get coordinates to cell in staggered grid
		xIdx = getBinIdx(xPoints, cellLen, point(0));
		yIdx = getBinIdx(yPoints, cellLen, point(1));
		zIdx = getBinIdx(zPoints, cellLen, point(2));

		// TODO: interpolate
		qdotCol(i) = 0;
	}
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

