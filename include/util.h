#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Core>

void addToCol(Eigen::MatrixXd& matrix, int columnIndex, double value);

void addRows(Eigen::MatrixXd& matrix, const Eigen::MatrixXd& rows);

#endif
