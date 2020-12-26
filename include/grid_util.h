#include <Eigen/Core>
#include <Eigen/Geometry>

#include "grid.h"

//
// Inputs:
//   V    #V by 3 matrix of vertex representing a box
// Outputs:
//   box  AlignedBox constructed from V
void createAlignedBox(const Eigen::MatrixXd& V, Eigen::AlignedBox3d& box);


//
// Scale, rotate, and translate the given mesh cloud
// - Based on the transformation between the vertices' bounding
//   box and the given destination bounding box
// - Translation is based on the BottomLeftFloor (Eigen::AlignedBox corner)
// - Rotations are at most 90 degrees. The input parameters do not allow definitions of larger rotations
//
// Inputs:
//   V               #V by 3 matrix of vertices
//   newOrientation  desired bounding box that transformed V should be in
//   keepSize        boolean indicating whether Vout should be the same size as V
void transformVertices(Eigen::MatrixXd& V, const Eigen::AlignedBox3d &newOrientation,  bool keepSize);

// keepSize = false by default
void transformVertices(Eigen::MatrixXd& V, const Eigen::AlignedBox3d &newOrientation);

//
// Aligns given mesh to the axis planes. Keeps the size the same
//
// Inputs:
//   V         #vertices by 3 matrix of vertices
// Outputs:
//   alignedV  #vertices by 3 matrix of vertices transformed to align to axis planes
void alignToAxis(const Eigen::MatrixXd& V);

void convertGridToPoints(Grid grid, Eigen::MatrixXd& points);

// For a point inside the grid, return the bin index that the point belongs to
// - Returns the smaller of the surrounding borders
int getBinIdx(std::vector<double> binBorders, double binSize, double currLocation);

// For a point in world-space on the border of the grid, return the index of the grid point
int getPointIdx(std::vector<double>& sortedPoints, double interval, double currPoint, double epsilon);
