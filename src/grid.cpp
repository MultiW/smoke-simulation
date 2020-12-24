#include "grid.h"
#include "grid_util.h"
#include "util.h"

#include <stdio.h>

Grid::Grid() {}

Grid::Grid(size_t d1, size_t d2, size_t d3):
	grid(d1, d2, d3)
{
}

void Grid::setWorldPoints(const Eigen::MatrixXd& points, double cellSize)
{
	// Grid's sorted x, y, z values
	// Note: Eigen Matrix is stored in column-major order
	colToSortedVector(points.col(0), this->_x);
	colToSortedVector(points.col(1), this->_y);
	colToSortedVector(points.col(2), this->_z);

	double epsilon = cellSize / 10.0;

	int xIdx, yIdx, zIdx;
	Eigen::RowVector3d currPoint;
	for (int i = 0; i < points.rows(); i++)
	{
		currPoint = points.row(i);

		xIdx = getPointIdx(this->_x, cellSize, currPoint(0), epsilon);
		yIdx = getPointIdx(this->_y, cellSize, currPoint(1), epsilon);
		zIdx = getPointIdx(this->_z, cellSize, currPoint(2), epsilon);

		if (xIdx < grid.size(0) && yIdx < grid.size(1) && zIdx < grid.size(2))
		{
			grid(xIdx, yIdx, zIdx).gridPoint = Eigen::Vector3i(xIdx, yIdx, zIdx);
			grid(xIdx, yIdx, zIdx).worldPoint = Eigen::Vector3d(currPoint(0), currPoint(1), currPoint(2));
		}
	}

	this->cellSize = cellSize;
	this->points = points;
}

void Grid::setPointValues(const Eigen::VectorXd& newValues)
{
	Eigen::Vector3d gridIdx;
	for (int i = 0; i < newValues.rows(); i++)
	{
		gridIdx = mapTo3d(i, this->size(0), this->size(1), this->size(2));
		this->grid(gridIdx(0), gridIdx(1), gridIdx(2)).value = newValues(i);
	}
}

double Grid::interpolatePoint(const Eigen::RowVector3d point)
{
	int xi, yi, zi; // cell in which point is located
	int x, y, z; // distance from corner of box to point

	// get coordinates to cell surrounding current point
	xi = getBinIdx(this->_x, this->cellSize, point(0));
	yi = getBinIdx(this->_y, this->cellSize, point(1));
	zi = getBinIdx(this->_z, this->cellSize, point(2));

	x = point(0) - this->grid(xi, yi, zi).worldPoint(0);
	y = point(1) - this->grid(xi, yi, zi).worldPoint(1);
	z = point(2) - this->grid(xi, yi, zi).worldPoint(2);

	// trilinear interpolation from all corners of current cell
	return this->safeGet(xi, yi, zi) * (this->cellSize - x) * (this->cellSize - y) * (this->cellSize - z)
		+ this->safeGet(xi + 1, yi, zi) * x * (this->cellSize - y) * (this->cellSize - z)
		+ this->safeGet(xi, yi + 1, zi) * (this->cellSize - x) * y * (this->cellSize - z)
		+ this->safeGet(xi, yi, zi + 1) * (this->cellSize - x) * (this->cellSize - y) * z
		+ this->safeGet(xi + 1, yi, zi + 1) * x * (this->cellSize - y) * z
		+ this->safeGet(xi, yi + 1, zi + 1) * (this->cellSize - x) * y * z
		+ this->safeGet(xi + 1, yi + 1, zi) * x * y * (this->cellSize - z)
		+ this->safeGet(xi + 1, yi + 1, zi + 1) * x * y * z;
}

void Grid::interpolatePoints(const Eigen::MatrixXd& q, Eigen::Ref<Eigen::VectorXd> qdotCol)
{
	for (int i = 0; i < q.rows(); i++)
	{
		qdotCol(i) = this->interpolatePoint(q.row(i));
	}
}

void Grid::setGridValues(const Eigen::MatrixXd& q, const Eigen::VectorXd qdotCol)
{
	this->clearValues();

	Eigen::RowVector3d point;
	int xi, yi, zi; // cell in which point is located
	int x, y, z; // distance from corner of box to point
	for (int i = 0; i < q.rows(); i++)
	{
		point = q.row(i);

		// get coordinates to cell surrounding current point
		xi = getBinIdx(this->_x, this->cellSize, point(0));
		yi = getBinIdx(this->_y, this->cellSize, point(1));
		zi = getBinIdx(this->_z, this->cellSize, point(2));

		x = point(0) - this->grid(xi, yi, zi).worldPoint(0);
		y = point(1) - this->grid(xi, yi, zi).worldPoint(1);
		z = point(2) - this->grid(xi, yi, zi).worldPoint(2);

		// update current cell's values with current point. Use weights from trilinear interpolation
		this->safeAdd(xi, yi, zi, qdotCol(i) * (this->cellSize - x) * (this->cellSize - y) * (this->cellSize - z));
		this->safeAdd(xi + 1, yi, zi, qdotCol(i) * x * (this->cellSize - y) * (this->cellSize - z));
		this->safeAdd(xi, yi + 1, zi, qdotCol(i) * (this->cellSize - x) * y * (this->cellSize - z));
		this->safeAdd(xi, yi, zi + 1, qdotCol(i) * (this->cellSize - x) * (this->cellSize - y) * z);
		this->safeAdd(xi + 1, yi, zi + 1, qdotCol(i) * x * (this->cellSize - y) * z);
		this->safeAdd(xi, yi + 1, zi + 1, qdotCol(i) * (this->cellSize - x) * y * z);
		this->safeAdd(xi + 1, yi + 1, zi, qdotCol(i) * x * y * (this->cellSize - z));
		this->safeAdd(xi + 1, yi + 1, zi + 1, qdotCol(i) * x * y * z);
	}
}


double Grid::safeGet(int i, int j, int k)
{
	if (i >= 0 && i < this->size(0) &&
		j >= 0 && j < this->size(1) &&
		k >= 0 && k < this->size(2))
	{
		return this->grid(i, j, k).value;
	}
	return 0;
}

void Grid::safeAdd(int i, int j, int k, double value)
{
	if (i >= 0 && i < this->size(0) &&
		j >= 0 && j < this->size(1) &&
		k >= 0 && k < this->size(2))
	{
		this->grid(i, j, k).value += value;
	}
}

void Grid::clearValues()
{
	this->setConstantValue(0);
}

void Grid::setConstantValue(double value)
{
	for (int i = 0; i < this->size(0); i++)
	{
		for (int j = 0; j < this->size(1); j++)
		{
			for (int k = 0; k < this->size(2); k++)
			{
				this->grid(i, j, k).value = value;
			}
		}
	}
}

std::vector<double> const& Grid::x()
{
	return this->_x;
}

std::vector<double> const& Grid::y()
{
	return this->_y;
}

std::vector<double> const& Grid::z() 
{
	return this->_z;
}

// ================================
// === List3d Wrapper Functions ===
// ================================

Point& Grid::operator()(size_t i, size_t j, size_t k)
{
	return this->grid(i, j, k);
}

Point const& Grid::operator()(size_t i, size_t j, size_t k) const
{
	return this->grid(i, j, k);
}


void Grid::resize(size_t d1, size_t d2, size_t d3)
{
	this->grid.resize(d1, d2, d3);
}

size_t Grid::size(int dimension)
{
	return this->grid.size(dimension);
}
