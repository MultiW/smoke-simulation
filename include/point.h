#include <Eigen/Core>

class Point {
public:
    Eigen::Vector3d worldPoint;
    Eigen::Vector3i gridPoint;
    double value;
};