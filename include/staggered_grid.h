#include <Eigen/Core>
#include "point.h"
#include "list3d.h"

class StaggeredGrid {

public:
	List3d<Point> uGrid;
	List3d<Point> vGrid;
	List3d<Point> wGrid;
	List3d<Point> pGrid;
	int dim;

	StaggeredGrid(int dim);
	
	void initVelocities(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot);

	void updateGrids();
};