#include "advection.h"

void advection(Eigen::VectorXd q, Eigen::VectorXd qdot, int dt) {
	q = q + qdot * dt;
}