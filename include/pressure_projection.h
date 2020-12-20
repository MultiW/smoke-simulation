#include <Eigen/Core>
#include <Eigen/Geometry>

void pressure_projection(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, Eigen::AlignedBox3d& boundaries, double pressure);