#ifndef GRID_H
#define GRID_H

#include "list3d.h"
#include "point.h"

#include <Eigen/Core>

#include <vector>

/*
* Wrapper for List3d<Point> class
*/
class Grid {
public:
	Grid();
	Grid::Grid(size_t d1, size_t d2, size_t d3, double cellSize);

	void setWorldPoints(const Eigen::MatrixXd& points);

	void setConstantValue(double value);

	void setRandomValues(double min, double max);

	void setYZPlane(int x, double value);
	void setXZPlane(int y, double value);
	void setXYPlane(int z, double value);

	/*
	* 3D grid of values flattened to a 1D array
	* Assumption: dimensions of given vector matches this grid's dimensions
	*/
	void setPointValues(const Eigen::VectorXd& newValues);

	/*
	* Get the trilinear interpolation value of the point from the enclosing cube's values.
	*/
	double interpolatePoint(const Eigen::RowVector3d point);
	double sharpInterpolatePoint(const Eigen::RowVector3d point);

	// sorted x, y, z values of the grid
	std::vector<double> const& x();
	std::vector<double> const& y();
	std::vector<double> const& z();

	// Returns the value of the identified point. Returns value of the nearest point if out of bounds
	double safeGet(int i, int j, int k);
	// Add to the value of the given point. Do nothing if given indices are out of bounds
	void safeAdd(int i, int j, int k, double value);

	// Boundary checks
	bool isInBounds(int i, int j, int k);
	bool isPointInBounds(double x, double y, double z);

	// Works if point is out of bounds
	void getNearestGridPoint(Eigen::Vector3i& gridPoint, const Eigen::RowVector3d& point);

	// Debugging
	void dumpValues();
	void dumpValues(int i, int j, int k);

	// Wrapper functions for List3d class
    Point& operator()(size_t i, size_t j, size_t k);
	Point const& operator()(size_t i, size_t j, size_t k) const;
	void resize(size_t d1, size_t d2, size_t d3);
	size_t size(int dimension);
private:
	Eigen::MatrixXd points; // raw points of the grid
	List3d<Point> grid;
	double cellSize;

	std::vector<double> _x;
	std::vector<double> _y;
	std::vector<double> _z;

	/* Set values of all points to 0 */
	void clearValues();
};

#endif