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

	// Smoke-specific components
	double ambientTemp;
	Grid tempGrid; // temperature
	Grid densityGrid;
public:
	StaggeredGrid();
	StaggeredGrid(const Eigen::AlignedBox3d& box, const Eigen::Vector3i& dim, double ambientTemp, double defaultDensity);

	void setGridVelocities(const Eigen::MatrixXd& q, const Eigen::MatrixXd& qdot);

	void updateSmokeAndDensity();

	void applyVorticityConfinement(const Eigen::MatrixXd& q, const Eigen::MatrixXd& qdot);

	/* 
	* 1. Compute pressure from the grid velocities
	* 2. Update grid velocities using pressure
	* 3. Interpolate particle velocities using grid velocities
	*/
	void computePressureProjections(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt);

	// For testing
	void getGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);
private:
	double getCellSize();

	void createGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);

	// Convert between particles velocities and grid velocities
	void getInterpolatedVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);

	void updateGridVelocities();
	void computePressure(Eigen::VectorXd p, double dt);
};