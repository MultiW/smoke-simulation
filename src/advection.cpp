#include "advection.h"

void advection(Eigen::MatrixXd &q, Eigen::MatrixXd &qdot, double dt) {
	q = q + qdot * dt;
}