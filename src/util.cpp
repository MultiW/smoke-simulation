#include "util.h"

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
