#include "grid.h"
#include "grid_util.h"
#include "util.h"

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
		else
		{
			// TODO: REMOVE. For debugging
			printf("Not adding point: %d %d %d\n", xIdx, yIdx, zIdx);
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

void Grid::interpolateToPoints(const Eigen::MatrixXd& q, Eigen::Ref<Eigen::VectorXd> qdotCol)
{
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

		// trilinear interpolation from all corners of current cell
		qdotCol(i) = 0; // clear velocity value
		qdotCol(i) += this->grid(xi, yi, zi).value * (this->cellSize - x) * (this->cellSize - y) * (this->cellSize - z);
		qdotCol(i) += this->grid(xi + 1, yi, zi).value * x * (this->cellSize - y) * (this->cellSize - z);
		qdotCol(i) += this->grid(xi, yi + 1, zi).value * (this->cellSize - x) * y * (this->cellSize - z);
		qdotCol(i) += this->grid(xi, yi, zi + 1).value * (this->cellSize - x) * (this->cellSize - y) * z;
		qdotCol(i) += this->grid(xi + 1, yi, zi + 1).value * x * (this->cellSize - y) * z;
		qdotCol(i) += this->grid(xi, yi + 1, zi + 1).value * (this->cellSize - x) * y * z;
		qdotCol(i) += this->grid(xi + 1, yi + 1, zi).value * x * y * (this->cellSize - z);
		qdotCol(i) += this->grid(xi + 1, yi + 1, zi + 1).value * x * y * z;
	}
}


void Grid::interpolateGridValues(const Eigen::MatrixXd& q, const Eigen::VectorXd qdotCol)
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
		this->grid(xi, yi, zi).value += qdotCol(i) * (this->cellSize - x) * (this->cellSize - y) * (this->cellSize - z);
		this->grid(xi + 1, yi, zi).value += qdotCol(i) * x * (this->cellSize - y) * (this->cellSize - z);
		this->grid(xi, yi + 1, zi).value += qdotCol(i) * (this->cellSize - x) * y * (this->cellSize - z);
		this->grid(xi, yi, zi + 1).value += qdotCol(i) * (this->cellSize - x) * (this->cellSize - y) * z;
		this->grid(xi + 1, yi, zi + 1).value += qdotCol(i) * x * (this->cellSize - y) * z;
		this->grid(xi, yi + 1, zi + 1).value += qdotCol(i) * (this->cellSize - x) * y * z;
		this->grid(xi + 1, yi + 1, zi).value += qdotCol(i) * x * y * (this->cellSize - z);
		this->grid(xi + 1, yi + 1, zi + 1).value += qdotCol(i) * x * y * z;
	}
}

void Grid::clearValues()
{
	for (int i = 0; i < this->size(0); i++)
	{
		for (int j = 0; j < this->size(1); j++)
		{
			for (int k = 0; k < this->size(2); k++)
			{
				this->grid(i, j, k).value = 0;
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
