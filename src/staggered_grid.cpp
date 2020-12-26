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
	const Eigen::MatrixXd& q,
	const Eigen::AlignedBox3d& box,
	const Eigen::Vector3i& dim,
	double cellSize
) :
	box(box),
	dim(dim),
	cellSize(cellSize),
	uGrid(Grid(dim(0), dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	vGrid(Grid(dim(0)- 1.0, dim(1), dim(2) - 1.0, cellSize)),
	wGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2), cellSize)),
	pGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),

	tempGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	densityGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),

	// vorticity confinement grids
	omegaUGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	omegaVGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	omegaWGrid(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	omegaNormal(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	omegaGradU(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	omegaGradV(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	omegaGradW(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	vortConfU(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	vortConfV(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize)),
	vortConfW(Grid(dim(0) - 1.0, dim(1) - 1.0, dim(2) - 1.0, cellSize))
{
	// Generate grids
	Eigen::MatrixXd u, v, w, cellCenters;
	this->createGridPoints(u, v, w, cellCenters);

	// Convert to Grid object
	this->uGrid.setWorldPoints(u);
	this->vGrid.setWorldPoints(v);
	this->wGrid.setWorldPoints(w);
	this->pGrid.setWorldPoints(cellCenters);

	this->tempGrid.setWorldPoints(cellCenters);
	this->densityGrid.setWorldPoints(cellCenters);

	this->omegaUGrid.setWorldPoints(cellCenters);
	this->omegaVGrid.setWorldPoints(cellCenters);
	this->omegaWGrid.setWorldPoints(cellCenters);
	this->omegaNormal.setWorldPoints(cellCenters);
	this->omegaGradU.setWorldPoints(cellCenters);
	this->omegaGradV.setWorldPoints(cellCenters);
	this->omegaGradW.setWorldPoints(cellCenters);
	this->vortConfU.setWorldPoints(cellCenters);
	this->vortConfV.setWorldPoints(cellCenters);
	this->vortConfW.setWorldPoints(cellCenters);

	this->initializeVelocities();
	this->initializeTemperatureAndDensity(q);
}

// ======================
// === Public classes ===
// ======================
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
	return this->cellSize;
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
	//this->uGrid.setRandomValues(0, this->getCellSize());
	//this->vGrid.setRandomValues(-this->getCellSize(), this->getCellSize());
	//this->wGrid.setRandomValues(-this->getCellSize(), this->getCellSize());
	double rand = 1;
	this->vGrid.setRandomValues(-rand, 0);
	this->uGrid.setRandomValues(-rand, 0);
	this->wGrid.setRandomValues(-rand, 0);

	// set velocities to not shoot outside the grid
	this->uGrid.setYZPlane(0, 0);
	this->uGrid.setYZPlane(this->uGrid.size(0) - 1, 0);
	this->vGrid.setXZPlane(0, 0);
	this->vGrid.setXZPlane(this->vGrid.size(1) - 1, 0);
	this->wGrid.setXYPlane(0, 0);
	this->wGrid.setXYPlane(this->wGrid.size(2) - 1, 0);
}


void StaggeredGrid::initializeTemperatureAndDensity(const Eigen::MatrixXd& q)
{
	for (int i = 0; i < q.rows(); i++)
	{
	}
	this->tempGrid.setConstantValue(1);
	this->densityGrid.setConstantValue(1);
}

// =================
// === Advection ===
// =================

// Given i, j, k of the current cell, return the velocity at the center of the grid 
// - Computed using faces of the surrounding cube
// - Faces of the cube that are out of bounds are defaulted to 0
void StaggeredGrid::getVelocityAtCenter(Eigen::Vector3d& velocityAtCenter, int i, int j, int k)
{
	velocityAtCenter(0) = (this->uGrid.safeGet(i, j, k) + this->uGrid.safeGet(i + 1, j, k)) / 2.0;
	velocityAtCenter(1) = (this->vGrid.safeGet(i, j, k) + this->vGrid.safeGet(i, j + 1, k)) / 2.0;
	velocityAtCenter(2) = (this->wGrid.safeGet(i, j, k) + this->wGrid.safeGet(i, j, k + 1)) / 2.0;
}

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
				this->getVelocityAtCenter(velocityAtCenter, i, j, k);
				prevPosition = grid(i, j, k).worldPoint - dt * velocityAtCenter;

				// The center position's temperature at the next timestep will be that imaginary particle's temperature
				// TODO: try removing this check, should be safe in theory
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
				// TODO: try removing this check, should be safe in theory
				if (grid.isPointInBounds(prevPos(0), prevPos(1), prevPos(2)))
				{
					grid(i, j, k).value = grid.interpolatePoint(prevPos);
				}
			}
		}
	}
}

void StaggeredGrid::advectPosition(Eigen::MatrixXd &q) {
	for (int i = 0; i < q.rows(); i++) {
		Eigen::RowVector3d point = q.row(i);
		Eigen::RowVector3d nextPoint, vel;
		this->getPointVelocity(vel, point);
		nextPoint = point + vel * dt;
		this->enforceBoundaries(point, nextPoint);
		q.row(i) = nextPoint;
	}
}

void StaggeredGrid::getPointVelocity(Eigen::RowVector3d &velocity, Eigen::RowVector3d &point) {
	velocity(0) = uGrid.interpolatePoint(point);
	velocity(1) = vGrid.interpolatePoint(point);
	velocity(2) = wGrid.interpolatePoint(point);
}

void StaggeredGrid::enforceBoundaries(const Eigen::RowVector3d &point, Eigen::RowVector3d &nextPoint) {
	Eigen::RowVector3d mins, maxs;
	mins << this->getMinX(), this->getMinY(), this->getMinZ();
	maxs << this->getMaxX(), this->getMaxY(), this->getMaxZ();

	Eigen::RowVector3d enclosedPoint = nextPoint;
	Eigen::RowVector3d currDist;
	double newDist, ratio;
	for (int i = 0; i < 3; i++) {
		if (enclosedPoint(i) < mins(i)) 
		{
			currDist = enclosedPoint - point;
			newDist = mins(i) - point(i); // new distance to boundary (that is less than the boundary)
			ratio = newDist / currDist(i); // ratio to trim distance to within boundary
			enclosedPoint = point + currDist * ratio;
		}
		else if (enclosedPoint(i) > maxs(i))
		{
			currDist = enclosedPoint - point;
			newDist = maxs(i) - point(i);
			ratio = newDist / currDist(i);
			enclosedPoint = point + currDist * ratio;
		}
	}

	nextPoint(0) = enclosedPoint(0);
	nextPoint(1) = enclosedPoint(1);
	nextPoint(2) = enclosedPoint(2);
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
	// - where alpha, beta are tuning constants`s
	// - s is the density and T is the temperature
	double fBuoy, s, T;	
	for (int i = 0; i < d1; i++)
	{
		for (int j = 0; j < d2; j++)
		{
			for (int k = 0; k < d3; k++)
			{
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
	double epsilon = 0.1;
	double dc = 0;
	this->vorticityConfinement(epsilon);

	for (int i = 1; i < uGrid.size(0); i++) {
		for (int j = 0; j < uGrid.size(1); j++) {
			for (int k = 0; k < uGrid.size(2); k++) {
				dc = dt * (this->vortConfU.safeGet(i, j, k) + this->vortConfU.safeGet(i + 1, j, k)) / (2 * FLUID_DENSITY);
				this->uGrid.safeAdd(i, j, k, dc);
			}
		}
	}


	for (int i = 1; i < vGrid.size(0); i++) {
		for (int j = 0; j < vGrid.size(1); j++) { 
			for (int k = 0; k < vGrid.size(2); k++) {
				dc = dt * (this->vortConfV.safeGet(i, j, k) + this->vortConfV.safeGet(i, j + 1, k)) / (2 * FLUID_DENSITY);
				this->vGrid.safeAdd(i, j, k, dc);
			}
		}
	}

	for (int i = 1; i < wGrid.size(0); i++) {
		for (int j = 0; j < wGrid.size(1); j++) {
			for (int k = 0; k < wGrid.size(2); k++) {
				dc = dt * (this->vortConfW.safeGet(i, j, k) + this->vortConfW.safeGet(i, j, k + 1)) / (2 * FLUID_DENSITY);
				this->wGrid.safeAdd(i, j, k, dc);
			}
		}
	}

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


void StaggeredGrid::centerVelandNorm(double denom) {
	for (int i = 0; i < this->omegaUGrid.size(0); i++) {
		for (int j = 0; j < this->omegaVGrid.size(1); j++) {
			for (int k = 0; k < this->omegaWGrid.size(2); k++) {

				//get the cetner velocities at each face and compute
				Eigen::Vector3d first, second, third, fourth;
				this->getVelocityAtCenter(first, i, j + 1, k);
				this->getVelocityAtCenter(second, i, j - 1, k);

				this->getVelocityAtCenter(third, i, j, k + 1);
				this->getVelocityAtCenter(fourth, i, j, k - 1);

				this->omegaUGrid(i, j, k).value = ((first[2] - second[2]) - (third[1] - fourth[1])) / denom;

				this->getVelocityAtCenter(first, i, j, k + 1);
				this->getVelocityAtCenter(second, i, j, k - 1);

				this->getVelocityAtCenter(third, i + 1, j, k);
				this->getVelocityAtCenter(fourth, i - 1, j, k);

				this->omegaVGrid(i, j, k).value = ((first[0] - second[0]) - (third[2] - fourth[2])) / denom;

				this->getVelocityAtCenter(first, i + 1, j, k);
				this->getVelocityAtCenter(second, i - 1, j, k);

				this->getVelocityAtCenter(third, i, j + 1, k);
				this->getVelocityAtCenter(fourth, i, j - 1, k);

				this->omegaWGrid(i, j, k).value = ((first[1] - second[1]) - (third[0] - fourth[0])) / denom;

				Eigen::Vector3d norm(this->omegaUGrid(i, j, k).value, this->omegaVGrid(i, j, k).value, this->omegaWGrid(i, j, k).value);
				this->omegaNormal(i, j, k).value = norm.norm();
			}

		}
	}
}

// Vorticity Confinement
void StaggeredGrid::vorticityConfinement(double epsilon) {
	double denom = 2 * this->getCellSize();

	this->centerVelandNorm(denom);

	for (int i = 0; i < this->omegaNormal.size(0); i++) {
		for (int j = 0; j < this->omegaNormal.size(1); j++) {
			for (int k = 0; k < this->omegaNormal.size(2); k++) {
				Eigen::Vector3d norm, omega, res;
				norm[0] = (this->omegaNormal.safeGet(i + 1, j, k) - this->omegaNormal.safeGet(i - 1, j, k)) / denom;
				norm[1] = (this->omegaNormal.safeGet(i, j + 1, k) - this->omegaNormal.safeGet(i, j - 1, k)) / denom;
				norm[2] = (this->omegaNormal.safeGet(i, j, k + 1) - this->omegaNormal.safeGet(i, j, k - 1)) / denom;
				norm = norm / (norm.norm() + 1e-20);
				omega << this->omegaUGrid.safeGet(i, j, k), this->omegaVGrid.safeGet(i, j, k), this->omegaWGrid.safeGet(i, j, k);

				res = norm.cross(omega) * this->getCellSize() * epsilon;
				
				this->vortConfU(i, j, k).value = res[0];
				this->vortConfV(i, j, k).value = res[1];
				this->vortConfW(i, j, k).value = res[2];

			}
		}
	}
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
