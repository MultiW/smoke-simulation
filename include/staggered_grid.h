#include <Eigen/Core>
#include <Eigen/Geometry>

#include "point.h"
#include "list3d.h"

typedef List3d<Point> Grid;

class StaggeredGrid {
	// Number of points along each dimension
	Eigen::Vector3i dim;

	// Dimension and location of grid in world-space
	Eigen::AlignedBox3d box;

	// Velocity grids
	Grid uGrid;
	Grid vGrid;
	Grid wGrid;
	
	// Pressure grid
	Grid pGrid;
public:
	StaggeredGrid();
	StaggeredGrid(const Eigen::AlignedBox3d& box, const Eigen::Vector3i& dim);
	void computeVelocity(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);

	// For testing
	void getGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);
private:
	double getCellSize();

	void updateWorldPoints(const Eigen::MatrixXd& points, Grid& grid);
	void createGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);

	void setGridVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);
	void getVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);
	void updateVelocityAndPressure();
};