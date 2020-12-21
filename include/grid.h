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

	void interpolateToPoints(const Eigen::MatrixXd& q, Eigen::Ref<Eigen::VectorXd> qdotCol);

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
};

#endif