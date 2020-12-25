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
	Grid::Grid(size_t d1, size_t d2, size_t d3);

	void setWorldPoints(const Eigen::MatrixXd& points, double cellSize);

	void setConstantValue(double value);

	/*
	* 3D grid of values flattened to a 1D array
	* Assumption: dimensions of given vector matches this grid's dimensions
	*/
	void setPointValues(const Eigen::VectorXd& newValues);

	/*
	* Get the trilinear interpolation value of the point from the enclosing cube's values.
	*/
	double interpolatePoint(const Eigen::RowVector3d point);

	/*
	* Set the grid values based on the weighted sum of the given points' values.
	* - Use trilinear interpolation weights for each point
	* Set to zero if no points are inside a cube
	*/
	void setGridValues(const Eigen::MatrixXd& q, const Eigen::VectorXd qdotCol);

	// sorted x, y, z values of the grid
	std::vector<double> const& x();
	std::vector<double> const& y();
	std::vector<double> const& z();

	/* Returns the value of the identified point. Return 0 if out of bounds */
	double safeGet(int i, int j, int k);

	/* Add to the value of the given point. Do nothing if given indices are out of bounds */
	void safeAdd(int i, int j, int k, double value);

	bool isInBounds(int i, int j, int k);
	bool isPointInBounds(double x, double y, double z);

	// Debugging
	void dumpValues();

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