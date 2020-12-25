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
	pGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0)),
	tempGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0)),
	densityGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0))
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

	this->initializeVelocities();

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

void StaggeredGrid::advectVelocities()
{
	advectVelocity(this->uGrid);
	advectVelocity(this->vGrid);
	advectVelocity(this->wGrid);
}

void StaggeredGrid::updateTemperatureAndDensity()
{
	this->advectCenterValues(this->tempGrid);
	this->advectCenterValues(this->densityGrid);
}

void StaggeredGrid::applyExternalForces()
{
	this->applyBuoyancyForce();
	this->applyVorticityConfinement();
}

void StaggeredGrid::applyPressureProjections() {
	Eigen::VectorXd p;
	this->computePressure(p);
	this->pGrid.setPointValues(p);

	// Update grid velocities using pressure
	this->updateGridVelocities();
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

// ======================
// === Initialization ===
// ======================

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


void StaggeredGrid::initializeVelocities()
{
	this->uGrid.setConstantValue(0);
	this->vGrid.setConstantValue(getRand(0, 1));
	this->wGrid.setConstantValue(0);
}

// =================
// === Advection ===
// =================

void StaggeredGrid::advectCenterValues(Grid& grid)
{
	int d1 = grid.size(0);
	int d2 = grid.size(1);
	int d3 = grid.size(2);

	Eigen::Vector3d prevPosition;
	Eigen::Vector3d velocityAtCenter;
	for (int i = 0; i < d1; i++) 
	{
		for (int j = 0; j < d2; j++) 
		{
			for (int k = 0; k < d3; k++) 
			{
				// Get location of the (imaginary) particle that will reach this center position after timestep dt
				velocityAtCenter(0) = (this->uGrid(i, j, k).value + this->uGrid(i + 1, j, k).value) / 2.0;
				velocityAtCenter(1) = (this->vGrid(i, j, k).value + this->vGrid(i, j + 1, k).value) / 2.0;
				velocityAtCenter(2) = (this->wGrid(i, j, k).value + this->wGrid(i, j, k + 1).value) / 2.0;
				prevPosition = grid(i, j, k).worldPoint - dt * velocityAtCenter;

				// The center position's temperature at the next timestep will be that imaginary particle's temperature
				if (grid.isPointInBounds(prevPosition(0), prevPosition(1), prevPosition(2)))
				{
					grid(i, j, k).value = grid.interpolatePoint(prevPosition);
				}
			}
		}
	}
}

void StaggeredGrid::advectVelocity(Grid& grid)
{
	int d1 = grid.size(0);
	int d2 = grid.size(1);
	int d3 = grid.size(2);

	Eigen::RowVector3d currPos;
	Eigen::RowVector3d currVel;
	Eigen::RowVector3d prevPos;
	for (int i = 0; i < d1; i++)
	{
		for (int j = 0; j < d2; j++)
		{
			for (int k = 0; k < d3; k++)
			{
				currPos = grid(i, j, k).worldPoint.transpose();
				this->getPointVelocity(currVel, currPos);
				prevPos = currPos - dt * currVel;
				grid(i, j, k).value = grid.interpolatePoint(prevPos);
			}
		}
	}
}

// =============================
// === Apply external forces ===
// =============================
void StaggeredGrid::applyBuoyancyForce()
{
	int d1 = this->vGrid.size(0);
	int d2 = this->vGrid.size(1);
	int d3 = this->vGrid.size(2);

	// Buoyancy force = (0, -alpha * s + beta(T - Tamb), 0)
	// - where alpha, beta are tuning constants
	// - s is the density and T is the temperature
	double fBuoy, s, T;	
	for (int i = 0; i < d1; i++)
	{
		for (int j = 0; j < d2; j++)
		{
			for (int k = 0; k < d3; k++)
			{
				// TODO: how to handle boundary cases? Only use one neighbor in those cases?
				//    - currently, the "out of bounds" neighbor is converted to a 0, 
				//      so the s and T values are small because we still divide by 2
				s = (this->densityGrid.safeGet(i, j, k) + this->densityGrid.safeGet(i, j - 1, k)) / 2;
				T = (this->tempGrid.safeGet(i, j, k) + this->tempGrid.safeGet(i, j - 1, k)) / 2;
				fBuoy = -ALPHA * s + BETA * (T - AMBIENT_TEMP);
				this->vGrid(i, j, k).value += dt * fBuoy;
			}
		}
	}
}

void StaggeredGrid::applyVorticityConfinement()
{
	// TODO: need more grids for vorticity, center velocity, etc.
	// Compute omega = curl of velocity
	Eigen::Vector3d omega;
}



// ======================================
// === Computing pressure projections ===
// ======================================
void StaggeredGrid::updateGridVelocities()
{
	double dp; // change in pressure
	double cellSize = this->getCellSize();

	int d1 = this->uGrid.size(0);
	int d2 = this->uGrid.size(1);
	int d3 = this->uGrid.size(2);
	for (int i = 1; i < d1 - 1; i++) {
		for (int j = 0; j < d2; j++) {
			for (int k = 0; k < d3; k++) {
				dp = (pGrid(i, j, k).value - pGrid(i - 1, j, k).value) / cellSize;
				uGrid(i, j, k).value -= dt / AIR_DENSITY * dp;
			}
		}
	}

	d1 = this->vGrid.size(0);
	d2 = this->vGrid.size(1);
	d3 = this->vGrid.size(2);

	for (int i = 0; i < d1; i++) {
		for (int j = 1; j < d2 - 1; j++) {
			for (int k = 0; k < d3; k++) {
				dp = (pGrid(i, j, k).value - pGrid(i, j - 1, k).value) / cellSize;
				vGrid(i, j, k).value -= dt / AIR_DENSITY * dp;
			}
		}
	}

	d1 = this->wGrid.size(0);
	d2 = this->wGrid.size(1);
	d3 = this->wGrid.size(2);

	for (int i = 0; i < d1; i++) {
		for (int j = 0; j < d2; j++) {
			for (int k = 1; k < d3 - 1; k++) {
				dp = (pGrid(i, j, k).value - pGrid(i, j, k - 1).value) / cellSize;
				wGrid(i, j, k).value -= dt / AIR_DENSITY * dp;
			}
		}
	}
}

void safe_push_back(std::vector<T>& vector, int row, int i, int j, int k, int d1, int d2, int d3, double val) {
	if (val != 0) {
		vector.push_back(T(row, mapTo1d(i, j, k, d1, d2, d3), val));
	}
}

void StaggeredGrid::computePressure(Eigen::VectorXd p)
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

// =================================
// === Simulation box boundaries ===
// =================================
double StaggeredGrid::getMinX()
{
	return this->box.corner(this->box.BottomLeftFloor).x();
}

double StaggeredGrid::getMinY()
{
	return this->box.corner(this->box.BottomLeftFloor).y();
}

double StaggeredGrid::getMinZ()
{
	return this->box.corner(this->box.BottomLeftFloor).z();
}

double StaggeredGrid::getMaxX()
{
	return this->box.corner(this->box.TopRightCeil).x();
}

double StaggeredGrid::getMaxY()
{
	return this->box.corner(this->box.TopRightCeil).y();
}

double StaggeredGrid::getMaxZ()
{
	return this->box.corner(this->box.TopRightCeil).z();
}
void StaggeredGrid::advectPosition(Eigen::MatrixXd &q) {
	//Eigen::MatrixXd qnext;
	//qnext.resize(q.rows(), q.cols());



	for (int i = 0; i < q.rows(); i++) {
		Eigen::RowVector3d point = q.row(i);
		Eigen::RowVector3d vel;
		this->getPointVelocity(vel, point);

		vel = vel * dt;
		this->enforceBoundaries(vel, point);

		//qnext.row(i) = point + vel*dt;

	}
}

void StaggeredGrid::getPointVelocity(Eigen::RowVector3d &velocity, Eigen::RowVector3d &point) {
	velocity(0) = uGrid.interpolatePoint(point);
	velocity(1) = vGrid.interpolatePoint(point);
	velocity(2) = wGrid.interpolatePoint(point);
}

void StaggeredGrid::enforceBoundaries(Eigen::RowVector3d &vel, Eigen::RowVector3d &point) {



	Eigen::RowVector3d mins;
	mins << this->getMinX(), this->getMinY(), this->getMinZ();
	Eigen::RowVector3d maxs;
	maxs << this->getMaxX(), this->getMaxY(), this->getMaxZ();
	Eigen::RowVector3d newPoint = vel + point;
	for (int i = 0; i < 3; i++) {
		if (newPoint[0] < mins[0]) {
			double dist = mins[0] - point[0];
			vel = vel / vel[i] * dist;
		}
		else if (newPoint[0] > maxs[1]) {
			double dist = maxs[0] - point[0];
			vel = vel / vel[i] * dist;
		}
	}
}