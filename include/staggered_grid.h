#include <Eigen/Core>
#include <Eigen/Geometry>

#include "point.h"
#include "list3d.h"

typedef List3d<Point> Grid;

class StaggeredGrid {
	// Number of points along each dimension
	Eigen::Vector3i dim;

	// Velocity grids
	Grid uGrid;
	Grid vGrid;
	Grid wGrid;
	
	// Pressure grid
	Grid pGrid;
public:
	StaggeredGrid();
	StaggeredGrid(Eigen::Vector3i dim);
	void computeVelocity(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);
private:
	void setGridVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);
	void getVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);
	void updateVelocityAndPressure();
};