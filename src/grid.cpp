#include "grid.h"
#include "util.h"

#include <vector>

Grid::Grid() {}

Grid::Grid(size_t d1, size_t d2, size_t d3):
	grid(d1, d2, d3)
{
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

void Grid::setWorldPoints(const Eigen::MatrixXd& points, double cellSize)
{
	// Grid's sorted x, y, z values
	// Note: Eigen Matrix is stored in column-major order
	std::vector<double> xPoints;
	std::vector<double> yPoints;
	std::vector<double> zPoints;
	colToSortedVector(points.col(0), xPoints);
	colToSortedVector(points.col(1), yPoints);
	colToSortedVector(points.col(2), zPoints);

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
