#include <Eigen/Core>
#include <Eigen/Geometry>

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
//   V                #V by 3 matrix of vertices
//   currOrientation  bounding box of V
//   newOrientation   desired bounding box that transformed V should be in
//   keepSize         boolean indicating whether Vout should be the same size as V
void transformVertices(Eigen::MatrixXd& V, const Eigen::AlignedBox3d& currOrientation, const Eigen::AlignedBox3d& newOrientation, bool keepSize);

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
