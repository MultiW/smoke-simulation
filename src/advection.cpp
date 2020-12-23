#include "advection.h"

void advection(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, double dt) {
	//TODO implement boundaries if necessary
	q = q + qdot * dt;
}