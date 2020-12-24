#include <Eigen/Core>
#include <Eigen/Geometry>

#include "grid.h"

class StaggeredGrid {
	// Number of points along each dimension
	Eigen::Vector3i dim;

	// Dimension and location of grid in world-space
	Eigen::AlignedBox3d box;
public:
	// TODO: move back to private when done
	// Velocity grids
	Grid uGrid;
	Grid vGrid;
	Grid wGrid;
	
	// Pressure grid
	Grid pGrid;

	// Smoke-specific components
	Grid tempGrid; // temperature
	Grid densityGrid; // fluid density

	StaggeredGrid();
	StaggeredGrid(const Eigen::AlignedBox3d& box, const Eigen::Vector3i& dim);

	void setGridVelocities(const Eigen::MatrixXd& q, const Eigen::MatrixXd& qdot);

	void advectVelocities();

	/* Update the temperature and density fields using the velocity field (velocity grids) */
	void updateTemperatureAndDensity();

	void applyExternalForces();

	/* 
	* 1. Compute pressure from the grid velocities
	* 2. Update grid velocities using pressure
	*/
	void applyPressureProjections();

	// For testing
	void getGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);
private:
	double getCellSize();

	void createGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);
	void initializeVelocities();

	// External forces
	void applyBuoyancyForce();
	void applyVorticityConfinement();

	/* Convert between particles velocities and grid velocities */
	void getInterpolatedVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot);

	/* Update the value (e.g. temperature) at the center cell based on the velocity field */
	void advectCenterValues(Grid& grid);

	void updateGridVelocities();

	void computePressure(Eigen::VectorXd p);

	// Boundaries to simulation
	double getMinX();
	double getMinY();
	double getMinZ();
	double getMaxX();
	double getMaxY();
	double getMaxZ();
};