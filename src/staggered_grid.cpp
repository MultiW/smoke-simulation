#include "staggered_grid.h"
#include "util.h"
#include "grid_util.h"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>


#include <igl/grid.h>

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <set>

typedef Eigen::Triplet<double> T;

StaggeredGrid::StaggeredGrid() {}

StaggeredGrid::StaggeredGrid(const Eigen::AlignedBox3d& box, const Eigen::Vector3i& dim) :
	box(box),
	dim(dim),
	uGrid(Grid(dim(0), dim(1) - 1.0, dim(2) - 1.0)),
	vGrid(Grid(dim(0)- 1.0, dim(1), dim(2) - 1.0)),
	wGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2))),
	pGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0))
{
	// World-space locations of grid-points
	Eigen::MatrixXd u, v, w, p;
	this->createGridPoints(u, v, w, p);

	// Update world points of grids
	this->uGrid.setWorldPoints(u, this->getCellSize());
	this->vGrid.setWorldPoints(v, this->getCellSize());
	this->wGrid.setWorldPoints(w, this->getCellSize());
	this->pGrid.setWorldPoints(p, this->getCellSize());
}

double StaggeredGrid::getCellSize()
{
	return this->box.sizes()(0) / (this->dim(0) - 1.0);
}

// =====================================
// === Setting world point locations ===
// =====================================

void StaggeredGrid::createGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p)
{
	Eigen::MatrixXd model;
	igl::grid(this->dim, model);
	transformVertices(model, this->box);

	double cellHalfLen = this->getCellSize() / 2.0;

	u = model;
	transformVertices(u, this->box);
	addToCol(u, 1, cellHalfLen);
	addToCol(u, 2, cellHalfLen);

	v = model;
	transformVertices(v, this->box);
	addToCol(v, 0, cellHalfLen);
	addToCol(v, 2, cellHalfLen);

	w = model;
	transformVertices(w, this->box);
	addToCol(w, 0, cellHalfLen);
	addToCol(w, 1, cellHalfLen);

	p = model;
	transformVertices(p, this->box);
	addToCol(p, 0, cellHalfLen);
	addToCol(p, 1, cellHalfLen);
	addToCol(p, 2, cellHalfLen);
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

void StaggeredGrid::getGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p)
{
	convertGridToPoints(this->uGrid, u);
	convertGridToPoints(this->vGrid, v);
	convertGridToPoints(this->wGrid, w);
	convertGridToPoints(this->pGrid, p);
}


// ======================================
// === Computing pressure projections ===
// ======================================

void StaggeredGrid::setGridVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) {
	//TODO interpolate qdot into grids
}

// binBorders - sorted array of boundaries of all bins
int getBinIdx(std::vector<double> binBorders, double binSize, double currLocation)
{
	int min = binBorders.front();
	int max = binBorders.back();
	assert(currLocation >= min && currLocation <= max);

	for (int i = 1; i < binBorders.size(); i++)
	{
		if (currLocation < binBorders[i])
		{
			return i - 1;
		}
	}

	printf("Error in getBinIdx(): could not find bin index.");
	throw;
}

void StaggeredGrid::interpolateGrid(Eigen::MatrixXd& q, Eigen::VectorXd& qdotCol, Grid& grid)
{
	double cellLen = this->getCellSize();

	// TODO: fix this. vector is wrong
	std::vector<double> xPoints;
	std::vector<double> yPoints;
	std::vector<double> zPoints;
	colToSortedVector(q.col(0), xPoints);
	colToSortedVector(q.col(1), yPoints);
	colToSortedVector(q.col(2), zPoints);

	Eigen::RowVector3d point;
	int xIdx, yIdx, zIdx; // cell in which point is located
	for (int i = 0; i < q.rows(); i++)
	{
		point = q.row(i);

		// get coordinates to cell in staggered grid
		xIdx = getBinIdx(xPoints, cellLen, point(0));
		yIdx = getBinIdx(yPoints, cellLen, point(1));
		zIdx = getBinIdx(zPoints, cellLen, point(2));

		// TODO: interpolate
		qdotCol(i) = 0;
	}
}

void StaggeredGrid::getVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) {
	// TODO: trillinear interpolation
}


void StaggeredGrid::updateVelocityAndPressure(double dt, double density) {
	// Update the pressure
	
	//update parameters where this is called for dt and density

	int gridSize = this->dim(0) * this->dim(1) * this->dim(2);

	double gridGranularity = getCellSize();
	Eigen::MatrixXd D;
	Eigen::VectorXd pj, B, qj;
	Eigen::VectorXd ;
	double invx = 1 / gridGranularity;
	double invy = 1 / gridGranularity;
	double invz = 1 / gridGranularity;
	B.resize(6);
	B << -invx, invx, -invy, invy, -invz, invz;
	pj.resize(7);
	D.resize(6, 7);
	D.setZero();
	D << -invx, 0, invx, 0, 0, 0, 0,
		0, invx, -invx, 0, 0, 0, 0,
		0, 0, invy, -invy, 0, 0, 0,
		0, 0, -invy, 0, invy, 0, 0,
		0, 0, invz, 0, 0, -invz, 0,
		0, 0, -invz, 0, 0, 0, invz;
	Eigen::VectorXd q, Aj,p;

	Aj.resize(7);
	Aj.setZero();
	Aj = B.transpose()* D;

	Eigen::SparseMatrix<double> A;
	std::vector<T> coefficients;

	A.resize(gridSize, gridSize);
	A.setZero();

	q.resize(gridSize);
	q.setZero();
	qj.resize(6);
	qj.setZero();

	for (int i = 0; i < dim(0) - 1; i++) {
		for (int j = 0; j < dim(1) - 1; j++) {
			for (int k = 0; k < dim(2) - 1; k++) {
				qj << uGrid(i, j, k).value, uGrid(i + 1, j, k).value, vGrid(i, j, k).value, vGrid(i, j + 1, k).value, wGrid(i, j, k).value, wGrid(i, j, k + 1).value;
				//pj << pGrid(i - 1, j, k).value, pGrid(i + 1, j, k).value, pGrid(i, j, k).value, pGrid(i, j - 1, k).value, pGrid(i, j + 1, k).value, pGrid(i, j, k - 1).value, pGrid(i, j, k + 1).value;
				q(mapTo1d(i, j, k, this->dim(0), this->dim(1), this->dim(2))) = (B.transpose()*qj)(0)*density/dt;
				int Arow = mapTo1d(i, j, k, this->dim(0),this->dim(1), this->dim(2));
				//Figure out borders
				coefficients.push_back(T(Arow, mapTo1d(i - 1, j, k, this->dim(0), this->dim(1), this->dim(2)), Aj(0)));
				coefficients.push_back(T(Arow, mapTo1d(i + 1, j, k, this->dim(0), this->dim(1), this->dim(2)), Aj(1)));
				coefficients.push_back(T(Arow, mapTo1d(i, j, k, this->dim(0), this->dim(1), this->dim(2)), Aj(2)));
				coefficients.push_back(T(Arow, mapTo1d(i, j - 1, k, this->dim(0), this->dim(1), this->dim(2)), Aj(3)));
				coefficients.push_back(T(Arow, mapTo1d(i, j + 1, k, this->dim(0), this->dim(1), this->dim(2)), Aj(4)));
				coefficients.push_back(T(Arow, mapTo1d(i, j, k - 1, this->dim(0), this->dim(1), this->dim(2)), Aj(5)));
				coefficients.push_back(T(Arow, mapTo1d(i, j, k + 1, this->dim(0), this->dim(1), this->dim(2)), Aj(6)));
			}
		}
	}
	
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.compute(A);
	p.resize(gridSize);
	p.setZero();
	p = cg.solve(q);	
}

void StaggeredGrid::computeVelocity(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt, double density) {
	this->setGridVelocities(q, qdot);
	this->updateVelocityAndPressure(dt, density);
	this->getVelocities(q, qdot);
}

