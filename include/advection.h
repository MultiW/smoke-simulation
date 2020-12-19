//includes here
#include <Eigen/Core>
#include <Eigen/Geometry>

//input
//qdot - velocity (nx1)
//q - position (nx1)
//dt - time change
//output
//q updated position

void advection(Eigen::VectorXd q, Eigen::VectorXd qdot, int dt);