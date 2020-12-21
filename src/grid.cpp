#include "grid.h"
#include "util.h"

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

void Grid::interpolateToPoints(const Eigen::MatrixXd& q, Eigen::Ref<Eigen::VectorXd> qdotCol)
{
	Eigen::RowVector3d point;
	int xi, yi, zi; // cell in which point is located
	int x, y, z; // distance from corner of box to point
	for (int i = 0; i < q.rows(); i++)
	{
		point = q.row(i);

		// get coordinates to cell in staggered grid
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
