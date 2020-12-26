#include "grid.h"
#include "grid_util.h"
#include "util.h"

#include <stdio.h>
#include <iostream>

Grid::Grid() {}

Grid::Grid(size_t d1, size_t d2, size_t d3, double cellSize):
	grid(d1, d2, d3),
	cellSize(cellSize)
{
}

void Grid::setWorldPoints(const Eigen::MatrixXd& points)
{
	// Grid's sorted x, y, z values
	// Note: Eigen Matrix is stored in column-major order
	colToSortedVector(points.col(0), this->_x);
	colToSortedVector(points.col(1), this->_y);
	colToSortedVector(points.col(2), this->_z);

	double epsilon = this->cellSize / 10.0;

	int xIdx, yIdx, zIdx;
	Eigen::RowVector3d currPoint;
	for (int i = 0; i < points.rows(); i++)
	{
		currPoint = points.row(i);

		xIdx = getPointIdx(this->_x, this->cellSize, currPoint(0), epsilon);
		yIdx = getPointIdx(this->_y, this->cellSize, currPoint(1), epsilon);
		zIdx = getPointIdx(this->_z, this->cellSize, currPoint(2), epsilon);

		if (xIdx < grid.size(0) && yIdx < grid.size(1) && zIdx < grid.size(2))
		{
			grid(xIdx, yIdx, zIdx).gridPoint = Eigen::Vector3i(xIdx, yIdx, zIdx);
			grid(xIdx, yIdx, zIdx).worldPoint = Eigen::Vector3d(currPoint(0), currPoint(1), currPoint(2));
		}
	}
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

	x = (point(0) - this->grid(xi, yi, zi).worldPoint(0)) / this->cellSize;
	y = (point(1) - this->grid(xi, yi, zi).worldPoint(1)) / this->cellSize;
	z = (point(2) - this->grid(xi, yi, zi).worldPoint(2)) / this->cellSize;

	// trilinear interpolation from all corners of current cell
	return this->safeGet(xi, yi, zi) * (1 - x) * (1 - y) * (1 - z)
		+ this->safeGet(xi + 1, yi, zi) * x * (1 - y) * (1 - z)
		+ this->safeGet(xi, yi + 1, zi) * (1 - x) * y * (1 - z)
		+ this->safeGet(xi, yi, zi + 1) * (1 - x) * (1 - y) * z
		+ this->safeGet(xi + 1, yi, zi + 1) * x * (1 - y) * z
		+ this->safeGet(xi, yi + 1, zi + 1) * (1 - x) * y * z
		+ this->safeGet(xi + 1, yi + 1, zi) * x * y * (1 - z)
		+ this->safeGet(xi + 1, yi + 1, zi + 1) * x * y * z;
}

double Grid::safeGet(int i, int j, int k)
{
	// clamp i, j, k to the range.
	i = std::min(std::max(i, 0), (int) this->size(0) - 1);
	j = std::min(std::max(j, 0), (int) this->size(1) - 1);
	k = std::min(std::max(k, 0), (int) this->size(2) - 1);
	return this->grid(i, j, k).value;
}

void Grid::safeAdd(int i, int j, int k, double value)
{
	if (this->isInBounds(i, j, k))
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


void Grid::setRandomValues(double min, double max)
{
	for (int i = 0; i < this->size(0); i++)
	{
		for (int j = 0; j < this->size(1); j++)
		{
			for (int k = 0; k < this->size(2); k++)
			{
				this->grid(i, j, k).value = getRand(min, max);;
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


// ============================
// === Boundary and Indexes ===
// ============================
void Grid::getNearestGridPoint(Eigen::Vector3i& gridPoint, const Eigen::RowVector3d& point)
{
	gridPoint(0) = getBinIdx(this->_x, this->cellSize, point(0));
	gridPoint(1) = getBinIdx(this->_y, this->cellSize, point(1));
	gridPoint(2) = getBinIdx(this->_z, this->cellSize, point(2));
}

bool Grid::isInBounds(int i, int j, int k)
{
	return i >= 0 && i < this->size(0) &&
		j >= 0 && j < this->size(1) &&
		k >= 0 && k < this->size(2);
}


bool Grid::isPointInBounds(double x, double y, double z)
{
	return x >= this->_x.front() && x <= this->_x.back() &&
		y >= this->_y.front() && y <= this->_y.back() &&
		z >= this->_z.front() && z <= this->_z.back();
}

// =================
// === Debugging ===
// =================
void Grid::dumpValues()
{
	printf("\n");
	for (int i = 0; i < this->size(0); i++)
	{
		for (int j = 0; j < this->size(1); j++)
		{
			for (int k = 0; k < this->size(2); k++)
			{
				printf("%f  ", this->grid(i, j, k).value);
			}
		}
	}
	printf("\n");
}

// ====================================
// === Set Values in Selected Plane ===
// ====================================

void Grid::setYZPlane(int x, double value)
{
	for (int j = 0; j < this->size(1); j++)
	{
		for (int k = 0; k < this->size(2); k++)
		{
			this->grid(x, j, k).value = value;
		}
	}
}

void Grid::setXZPlane(int y, double value)
{
	for (int i = 0; i < this->size(0); i++)
	{
		for (int k = 0; k < this->size(2); k++)
		{
			this->grid(i, y, k).value = value;
		}
	}

}

void Grid::setXYPlane(int z, double value)
{
	for (int i = 0; i < this->size(0); i++)
	{
		for (int j = 0; j < this->size(1); j++)
		{
			this->grid(i, j, z).value = value;
		}
	}
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
