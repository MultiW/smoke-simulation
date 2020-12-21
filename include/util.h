#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Core>

#include <vector>

void addToCol(Eigen::MatrixXd& matrix, int columnIndex, double value);

void addRows(Eigen::MatrixXd& matrix, const Eigen::MatrixXd& rows);

void colToSortedVector(const Eigen::VectorXd& column, std::vector<double>& sortedVector);

#endif
