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

// ======================
// === Public classes ===
// ======================

void StaggeredGrid::updateSimulation(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt, double density) {
	// TODO: uncomment when all have been implemented
	//this->setGridVelocities(q, qdot);
	//this->updateVelocityAndPressure(dt, density);
	//this->getInterpolatedVelocities(q, qdot);
}

double StaggeredGrid::getCellSize()
{
	return this->box.sizes()(0) / (this->dim(0) - 1.0);
}

void StaggeredGrid::getGridPoints(Eigen::MatrixXd& u, Eigen::MatrixXd& v, Eigen::MatrixXd& w, Eigen::MatrixXd& p)
{
	convertGridToPoints(this->uGrid, u);
	convertGridToPoints(this->vGrid, v);
	convertGridToPoints(this->wGrid, w);
	convertGridToPoints(this->pGrid, p);
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

// ======================================
// === Computing pressure projections ===
// ======================================

void StaggeredGrid::setGridVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) 
{
	this->uGrid.interpolateGridValues(q, qdot.col(0));
	this->vGrid.interpolateGridValues(q, qdot.col(1));
	this->wGrid.interpolateGridValues(q, qdot.col(2));
}

void StaggeredGrid::getInterpolatedVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) 
{
	this->uGrid.interpolateToPoints(q, qdot.col(0));
	this->vGrid.interpolateToPoints(q, qdot.col(1));
	this->wGrid.interpolateToPoints(q, qdot.col(2));
}

void StaggeredGrid::updateVelocityAndPressure(double dt, double density) 
{
	Eigen::VectorXd p;
	this->computePressure(p, dt, density);

	// apply new pressure values to pGrid
	this->pGrid.setPointValues(p);

	// TODO: update velocity based on pressure
	this->updateGridVelocities();
}

void StaggeredGrid::updateGridVelocities()
{

}

void StaggeredGrid::computePressure(Eigen::VectorXd p, double dt, double density)
{
	// Number of cells in staggered grid (aka number of pressure points)
	int d1 = this->pGrid.size(0);
	int d2 = this->pGrid.size(1);
	int d3 = this->pGrid.size(2);
	int gridSize = d1 * d2 * d3;

	double gridGranularity = getCellSize();
	double invx = 1 / gridGranularity;
	double invy = 1 / gridGranularity;
	double invz = 1 / gridGranularity;

	// == Solve for pj (for every cell): B * D * pj = (density/dt) * B * qj ==
	//   - assemble globally to solve A * p = f

	// same for all cells: B, D, Aj = B * D
	Eigen::RowVectorXd B;
	B.resize(6);
	B << -invx, invx, -invy, invy, -invz, invz;

	Eigen::MatrixXd D;
	D.resize(6, 7);
	D.setZero();
	D << -invx, 0, invx, 0, 0, 0, 0,
		0, invx, -invx, 0, 0, 0, 0,
		0, 0, invy, -invy, 0, 0, 0,
		0, 0, -invy, 0, invy, 0, 0,
		0, 0, invz, 0, 0, -invz, 0,
		0, 0, -invz, 0, 0, 0, invz;

	Eigen::VectorXd Aj;
	Aj.resize(7);
	Aj.setZero();
	Aj = B * D;

	// global equation variable: p, f
	Eigen::VectorXd f;
	p.resize(gridSize);
	p.setZero();
	f.resize(gridSize);
	f.setZero();

	Eigen::SparseMatrix<double> A;
	A.resize(gridSize, gridSize);
	A.setZero();
	std::vector<T> coefficients;

	// TODO: remove pj
	// per-cell variables: pj, fj
	Eigen::VectorXd pj, fj;
	fj.resize(6);
	fj.setZero();
	pj.resize(7);

	for (int i = 0; i < d1; i++) {
		for (int j = 0; j < d2; j++) {
			for (int k = 0; k < d3; k++) {
				fj << 
					uGrid(i, j, k).value, 
					uGrid(i + 1, j, k).value, 
					vGrid(i, j, k).value, 
					vGrid(i, j + 1, k).value, 
					wGrid(i, j, k).value, 
					wGrid(i, j, k + 1).value;

				//pj << pGrid(i - 1, j, k).value, pGrid(i + 1, j, k).value, pGrid(i, j, k).value, pGrid(i, j - 1, k).value, pGrid(i, j + 1, k).value, pGrid(i, j, k - 1).value, pGrid(i, j, k + 1).value;

				// Assemble to global A, f 
				int row = mapTo1d(i, j, k, d1, d2, d3);
				f(row) = (B * fj)(0) * density / dt;

				//TODO: Figure out borders
				coefficients.push_back(T(row, mapTo1d(i - 1, j, k, d1, d2, d3), Aj(0)));
				coefficients.push_back(T(row, mapTo1d(i + 1, j, k, d1, d2, d3), Aj(1)));
				coefficients.push_back(T(row, mapTo1d(i, j, k, d1, d2, d3), Aj(2)));
				coefficients.push_back(T(row, mapTo1d(i, j - 1, k, d1, d2, d3), Aj(3)));
				coefficients.push_back(T(row, mapTo1d(i, j + 1, k, d1, d2, d3), Aj(4)));
				coefficients.push_back(T(row, mapTo1d(i, j, k - 1, d1, d2, d3), Aj(5)));
				coefficients.push_back(T(row, mapTo1d(i, j, k + 1, d1, d2, d3), Aj(6)));
			}
		}
	}
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	// Solve for pressure
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.compute(A);
	p = cg.solve(f);
}

