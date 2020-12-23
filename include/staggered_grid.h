#include <Eigen/Core>
#include <Eigen/Geometry>

#include "grid.h"

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

	/* 
	* 1. Update grid velocities with particle velocities
	* 2. Compute pressure using grid velocities
	* 3. Update grid velocities using pressure
	* 4. Interpolate particle velocities using grid velocities
	*/
	void computePressureProjections(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt, double density);

	void setGridVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);

	// For testing
	void getGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);
private:
	double getCellSize();

	void createGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);

	// Convert between particles velocities and grid velocities
	void getInterpolatedVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);

	void updateVelocityAndPressure(double dt, double density);
	void updateGridVelocities();
	void computePressure(Eigen::VectorXd p, double dt, double density);
};