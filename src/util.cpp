#include "util.h"

#include <set>
#include <algorithm>

void addToCol(Eigen::MatrixXd& matrix, int columnIndex, double value)
{
	Eigen::VectorXd ones;
	ones.resize(matrix.rows(), 1);
	ones.setOnes();
	matrix.col(columnIndex) += ones * value;
}

void addRows(Eigen::MatrixXd& matrix, const Eigen::MatrixXd& rows)
{
	assert(matrix.cols() == rows.cols());

	matrix.conservativeResize(matrix.rows() + rows.rows(), matrix.cols());
	matrix.block(matrix.rows() - rows.rows(), 0, rows.rows(), rows.cols()) = rows;
}

void colToSortedVector(const Eigen::VectorXd& column, std::vector<double>& sortedVector)
{
	std::set<double> uniqueValues(column.data(), column.data() + column.rows());
	sortedVector.assign(uniqueValues.begin(), uniqueValues.end());
	std::sort(sortedVector.begin(), sortedVector.end());
}

int mapTo1d(int i, int j, int k, int d1, int d2, int d3) {
	return i * d2 * d3 + j * d3 + k;
}

Eigen::Vector3d mapTo3d(int ind, int d1, int d2, int d3) {
	Eigen::Vector3d indices;

	indices(0) = ind / (d2 * d3);
	indices(1) = (ind % (d2 * d3)) / d3;
	indices(2) = ind - indices(0) * d2 * d3 - indices(1) * d3;
	assert(indices(2) == (ind % (d2 * d3)) % d3);
	return indices;
}

void flatten3d(Eigen::VectorXd& newArray, Grid& grid) {
	Eigen::Vector3i dim;
	dim << grid.size(0), grid.size(1), grid.size(2);
	newArray.resize(dim(0) * dim(1) * dim(2));
	newArray.setZero();
	for (int i = 0; i < dim(0) - 1; i++) {
		for (int j = 0; j < dim(1) - 1; j++) {
			for (int k = 0; k < dim(2) - 1; k++) {
				newArray(i * dim(1) * dim(2) + j * dim(2) + k) = grid(i, j, k).value;
			}
		}
	}
}

void unflatten(Grid& grid, Eigen::VectorXd& vector) {
	Eigen::Vector3i dim;
	dim << grid.size(0), grid.size(1), grid.size(2);
	for (int i = 0; i < dim(0) - 1; i++) {
		for (int j = 0; j < dim(1) - 1; j++) {
			for (int k = 0; k < dim(2) - 1; k++) {
				grid(i, j, k).value = vector(i * dim(1) * dim(2) + j * dim(2) + k);
			}
		}
	}
}
