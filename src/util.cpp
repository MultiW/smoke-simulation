#include "util.h"

inline void addToCol(Eigen::MatrixXd& matrix, int columnIndex, double value)
{
	Eigen::VectorXd ones;
	ones.resize(matrix.rows(), 1);
	ones.setOnes();
	matrix.col(columnIndex) += ones * value;
}
