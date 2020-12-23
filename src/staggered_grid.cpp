#include "staggered_grid.h"
#include "util.h"
#include "grid_util.h"
#include "constants.h"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <igl/grid.h>

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <set>

typedef Eigen::Triplet<double> T;


StaggeredGrid::StaggeredGrid() {}

StaggeredGrid::StaggeredGrid(
	const Eigen::AlignedBox3d& box,
	const Eigen::Vector3i& dim
) :
	box(box),
	dim(dim),
	uGrid(Grid(dim(0), dim(1) - 1.0, dim(2) - 1.0)),
	vGrid(Grid(dim(0)- 1.0, dim(1), dim(2) - 1.0)),
	wGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2))),
	pGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0))
{
	// Generate grids
	Eigen::MatrixXd u, v, w, cellCenters;
	this->createGridPoints(u, v, w, cellCenters);

	// Convert to Grid object
	this->uGrid.setWorldPoints(u, this->getCellSize());
	this->vGrid.setWorldPoints(v, this->getCellSize());
	this->wGrid.setWorldPoints(w, this->getCellSize());
	this->pGrid.setWorldPoints(cellCenters, this->getCellSize());
	this->tempGrid.setWorldPoints(cellCenters, this->getCellSize());
	this->densityGrid.setWorldPoints(cellCenters, this->getCellSize());

	// Set default temperature and density
	this->tempGrid.setConstantValue(FLUID_TEMP);
	this->densityGrid.setConstantValue(FLUID_DENSITY);
}

// ======================
// === Public classes ===
// ======================

void StaggeredGrid::setGridVelocities(const Eigen::MatrixXd& q, const Eigen::MatrixXd& qdot)
{
	this->uGrid.setGridValues(q, qdot.col(0));
	this->vGrid.setGridValues(q, qdot.col(1));
	this->wGrid.setGridValues(q, qdot.col(2));
}

void StaggeredGrid::updateSmokeAndDensity()
{
	// TODO: Xin
}

void StaggeredGrid::applyVorticityConfinement(const Eigen::MatrixXd& q, const Eigen::MatrixXd& qdot)
{
	// TODO: Xin
}

void StaggeredGrid::computePressureProjections(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot, double dt) {
	Eigen::VectorXd p;
	this->computePressure(p, dt);
	this->pGrid.setPointValues(p);

	// Update grid velocities using pressure
	this->updateGridVelocities();

	this->getInterpolatedVelocities(q, qdot);
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

// ========================================
// === Set grid and particle velocities ===
// ========================================

void StaggeredGrid::getInterpolatedVelocities(Eigen::MatrixXd& q, Eigen::MatrixXd& qdot) 
{
	this->uGrid.interpolatePoints(q, qdot.col(0));
	this->vGrid.interpolatePoints(q, qdot.col(1));
	this->wGrid.interpolatePoints(q, qdot.col(2));
}

// ======================================
// === Computing pressure projections ===
// ======================================

void StaggeredGrid::updateGridVelocities()
{
	// TODO: James

	int d1 = this->uGrid.size(0);
	int d2 = this->uGrid.size(1);
	int d3 = this->uGrid.size(2);

	for (int i = 0; i < d1; i++) {
		for (int j = 0; j < d2; j++) {
			for (int k = 0; k < d3; k++) {
				if (i == 0) {
					continue;
				}
				uGrid(i, j, k).value = pGrid(i, j, k).value - pGrid(i-1, j, k).value;
			}
		}
	}

	d1 = this->vGrid.size(0);
	d2 = this->vGrid.size(1);
	d3 = this->vGrid.size(2);

	for (int i = 0; i < d1; i++) {
		for (int j = 0; j < d2; j++) {
			for (int k = 0; k < d3; k++) {
				if (j == 0) {
					continue;
				}
				vGrid(i, j, k).value = pGrid(i, j, k).value - pGrid(i, j - 1, k).value;
			}
		}
	}

	d1 = this->wGrid.size(0);
	d2 = this->wGrid.size(1);
	d3 = this->wGrid.size(2);

	for (int i = 0; i < d1; i++) {
		for (int j = 0; j < d2; j++) {
			for (int k = 0; k < d3; k++) {
				if (k == 0) {
					continue;
				}
				wGrid(i, j, k).value = pGrid(i, j, k).value - pGrid(i, j, k - 1).value;
			}
		}
	}
}

void safe_push_back(std::vector<T>& vector, int row, int i, int j, int k, int d1, int d2, int d3, double val) {
	if (val != 0) {
		vector.push_back(T(row, mapTo1d(i, j, k, d1, d2, d3), val));
	}
}

void StaggeredGrid::computePressure(Eigen::VectorXd p, double dt)
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

	Eigen::MatrixXd D, selection;
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

	// per-cell variables: pj, fj
	Eigen::VectorXd qj;
	qj.resize(6);
	qj.setZero();

	selection.resize(6, 6);
	selection.setZero();
	selection.setIdentity();

	for (int i = 0; i < d1; i++) {
		for (int j = 0; j < d2; j++) {
			for (int k = 0; k < d3; k++) {
				qj << 
					uGrid(i, j, k).value, 
					uGrid(i + 1, j, k).value, 
					vGrid(i, j, k).value, 
					vGrid(i, j + 1, k).value, 
					wGrid(i, j, k).value, 
					wGrid(i, j, k + 1).value;

				selection.setIdentity();

				if (i == 0) {
					selection(0, 0) = 0;
				} else if(i == d1 - 1) {
					selection(1, 1) = 0;
				}

				if (j == 0) {
					selection(2, 2) = 0;
				}
				else if (j == d2 - 1) {
					selection(3, 3) = 0;
				}

				if (k == 0) {
					selection(4, 4) = 0;
				}
				else if (k == d3 - 1) {
					selection(5, 5) = 0;
				}

				Aj = B * selection * D;

				// Assemble to global A, f 
				int row = mapTo1d(i, j, k, d1, d2, d3);
				f(row) = (B * selection *qj)(0) * AIR_DENSITY / dt;

				safe_push_back(coefficients, row, i - 1, j, k, d1, d2, d3, Aj(0));
				safe_push_back(coefficients, row, i + 1, j, k, d1, d2, d3, Aj(1));
				safe_push_back(coefficients, row, i, j, k, d1, d2, d3, Aj(2));
				safe_push_back(coefficients, row, i, j - 1, k, d1, d2, d3, Aj(3));
				safe_push_back(coefficients, row, i, j + 1, k, d1, d2, d3, Aj(4));
				safe_push_back(coefficients, row, i, j, k - 1, d1, d2, d3, Aj(5));
				safe_push_back(coefficients, row, i, j, k + 1, d1, d2, d3, Aj(6));
			}
		}
	}

	A.setFromTriplets(coefficients.begin(), coefficients.end());

	// Solve for pressure
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.compute(A);
	p = cg.solve(f);
}


