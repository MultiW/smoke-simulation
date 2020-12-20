#include "grid_util.h"

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>
#include <igl/rotation_matrix_from_directions.h>

const static int GRID_LEN = 50; // number of voxels on each side of the grid

void createAlignedBox(const Eigen::MatrixXd &V, Eigen::AlignedBox3d &box) {
	for (int i = 0; i < V.rows(); i++) {
		box.extend(V.row(i).transpose());
	}
}

void transformVertices(Eigen::MatrixXd &V, const Eigen::AlignedBox3d &newOrientation, bool keepSize) {
	Eigen::Vector3d corner1 = newOrientation.corner(newOrientation.BottomLeftFloor);
	Eigen::Vector3d corner2 = newOrientation.corner(newOrientation.TopRightCeil);

	// Compute diagonals of input and output grids
	Eigen::AlignedBox<double, 3> inBox;
	createAlignedBox(V, inBox);
	Eigen::Vector3d inBLF = inBox.corner(inBox.BottomLeftFloor);
	Eigen::Vector3d inTRC = inBox.corner(inBox.TopRightCeil);
	Eigen::Vector3d inDiag = inBLF - inTRC;

	Eigen::Vector3d outDiag = corner1 - corner2;

	// Scale V to desired size
	Eigen::Matrix3d scale = Eigen::Matrix3d::Identity();
	if (!keepSize) {
		scale.diagonal() << outDiag.cwiseQuotient(inDiag).array();
	}
	V *= scale;

	// Rotate
	inDiag = scale * inBLF - scale * inTRC;
	Eigen::Matrix3d R = igl::rotation_matrix_from_directions(inDiag.normalized(), outDiag.normalized());
	V *= R.transpose();

	// Move V to desired location
	Eigen::Vector3d translate = corner1 - R * (scale * inBLF);
	for (int i = 0; i < V.rows(); i++) {
		V.row(i) += translate.transpose();
	}
}

void transformVertices(Eigen::MatrixXd& V, const Eigen::AlignedBox3d& newOrientation) {
	transformVertices(V, newOrientation, false);
}

void alignToAxis(Eigen::MatrixXd& V) {
	Eigen::AlignedBox3d inputBox;
	createAlignedBox(V, inputBox);

	Eigen::AlignedBox3d alignedBox;
	Eigen::MatrixXd alignedGrid;
	Eigen::Vector3d dim = Eigen::Vector3d::Constant(10.0);
	igl::grid(dim, alignedGrid);
	createAlignedBox(alignedGrid, alignedBox);

	transformVertices(V, alignedBox, true);
}