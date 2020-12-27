#include <Eigen/Core>
#include <Eigen/Geometry>

#include "grid.h"

class StaggeredGrid {
	// Number of points along each dimension
	Eigen::Vector3i dim;

	// Dimension and location of grid in world-space
	Eigen::AlignedBox3d box;

	double cellSize;

	// External objects
	Eigen::RowVector3d ballCenter;

	Eigen::MatrixXd* bunnyV;
	Eigen::MatrixXi* bunnyF;
public:
	// TODO: move back to private when done
	// Velocity grids
	Grid uGrid;
	Grid vGrid;
	Grid wGrid;
	
	// Pressure grid
	Grid pGrid;

	// Omega grid used in vorticity
	Grid omegaUGrid;
	Grid omegaVGrid;
	Grid omegaWGrid;

	Grid omegaNormal;
	Grid omegaGradU;
	Grid omegaGradV;
	Grid omegaGradW;

	Grid vortConfU;
	Grid vortConfV;
	Grid vortConfW;

	// Smoke-specific components
	Grid tempGrid; // temperature
	Grid densityGrid; // fluid density
	// -----------------------------------------

	StaggeredGrid();
	StaggeredGrid(const Eigen::MatrixXd& q, const Eigen::AlignedBox3d& box, const Eigen::Vector3i& dim, double cellSize);

	// Advection
	void advectVelocities();
	void advectPosition(Eigen::MatrixXd &q);

	// Update the temperature and density fields using the velocity field (velocity grids)
	void updateTemperatureAndDensity();

	void applyExternalForces();

	/* 
	* 1. Compute pressure from the grid velocities
	* 2. Update grid velocities using pressure
	*/
	void applyPressureProjections();

	//Vorticity Confinement
	void vorticityConfinement(double epislon);
	void centerVelandNorm(double denom);


	//for balls and bunnies
	void updateExternalObjects(Eigen::RowVector3d& ballCenter, Eigen::MatrixXd* bunnyV, Eigen::MatrixXi* bunnyF);

	// For testing
	void getGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);
private:
	double getCellSize();

	void getVelocityAtCenter(Eigen::Vector3d& velocityAtCenter, int i, int j, int k);

	// Initialization
	void createGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p);
	void initializeVelocities();
	void initializeTemperatureAndDensity(const Eigen::MatrixXd& q);

	// External forces
	void applyBuoyancyForce();
	void applyVorticityConfinement();

	// Advection
	void advectCenterValues(Grid& grid); // Update the value (e.g. temperature) at the center cell based on the velocity field
	void advectVelocity(Grid& grid);
	void enforceBoundaries(const Eigen::RowVector3d& point, Eigen::RowVector3d& nextPoint);
	void getPointVelocity(Eigen::RowVector3d &velocity, Eigen::RowVector3d &point);

	// Pressure projection
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