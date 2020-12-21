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

	/*
	* 3D grid of values flattened to a 1D array
	* Assumption: dimensions of given vector matches this grid's dimensions
	*/
	void Grid::setPointValues(const Eigen::VectorXd& newValues);

	/*
	* For each point q, interpolate its qdotCol values from enclosing cube's values
	*/
	void interpolateToPoints(const Eigen::MatrixXd& q, Eigen::Ref<Eigen::VectorXd> qdotCol);

	/*
	* Set the grid values based on the given points' values.
	* Set to zero if no points are inside a cube
	*/
	void interpolateGridValues(const Eigen::MatrixXd& q, const Eigen::VectorXd qdotCol);

	// sorted x, y, z values of the grid
	std::vector<double> const& x();
	std::vector<double> const& y();
	std::vector<double> const& z();

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