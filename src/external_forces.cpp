#include "external_forces.h"

// acceleration of gravity
const Eigen::RowVector3d gravity(0, -9.8, 0);

void external_forces(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt)
{
	for (int i = 0; i < qdot.rows(); i++)
	{
		qdot.row(i) += gravity * dt;
	}
}
