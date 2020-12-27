#include "grid_util.h"
#include "util.h"

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>
#include <igl/rotation_matrix_from_directions.h>

#include <stdio.h>

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

void convertGridToPoints(Grid grid, Eigen::MatrixXd& points)
{
	points.resize(0, 3);
	for (int x = 0; x < grid.size(0); x++)
	{
		for (int y = 0; y < grid.size(1); y++)
		{
			for (int z = 0; z < grid.size(2); z++)
			{
				addRows(points, grid(x, y, z).worldPoint.transpose());
			}
		}
	}
}

// binBorders - sorted array of boundaries of all bins
int getBinIdx(std::vector<double> binBorders, double binSize, double currLocation)
{
	double min = binBorders.front();
	double max = binBorders.back();

	if (currLocation <= min) {
		return 0;
	}
	else if (currLocation >= max) {
		return binBorders.size() - 1;
	}

	for (int i = 1; i < binBorders.size(); i++)
	{
		if (currLocation < binBorders[i])
		{
			return i - 1;
		}
	}
	printf("Error in getBinIdx(): could not find bin index.\n");
	throw;
}

int getPointIdx(std::vector<double> &sortedPoints, double interval, double currPoint, double epsilon)
{
	double min = sortedPoints.front();
	double max = sortedPoints.back();

	// Check that point is in bounds
	assert(currPoint >= min - epsilon && currPoint <= max + epsilon);

	int borderIdx = std::round((currPoint - min) / interval);

	// Check that point is within "epsilon" from nearest grid point
	assert(std::abs(currPoint - sortedPoints[borderIdx]) < epsilon);

	return borderIdx;
}

bool isInBox(const Eigen::AlignedBox3d& box, Eigen::Vector3d point)
{
	Eigen::Vector3d minCorner = box.corner(Eigen::AlignedBox3d::BottomLeftFloor);
	Eigen::Vector3d maxCorner = box.corner(Eigen::AlignedBox3d::TopRightCeil);
	return
		point(0) >= minCorner(0) && point(0) <= maxCorner(0) &&
		point(1) >= minCorner(1) && point(1) <= maxCorner(1) &&
		point(2) >= minCorner(2) && point(2) <= maxCorner(2);
}
