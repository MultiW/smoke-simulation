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
