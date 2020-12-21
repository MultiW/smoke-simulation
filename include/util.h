#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Core>

#include "grid.h"
#include <vector>

void addToCol(Eigen::MatrixXd& matrix, int columnIndex, double value);

void addRows(Eigen::MatrixXd& matrix, const Eigen::MatrixXd& rows);

void colToSortedVector(const Eigen::VectorXd& column, std::vector<double>& sortedVector);

void flatten3d(Eigen::VectorXd& newArray, Grid &grid);

void unflatten(Grid &grid, Eigen::VectorXd &vector);

#endif
