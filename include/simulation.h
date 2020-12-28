#ifndef SIMULATION_H
#define SIMULATION_H

#include "visualization.h"

#include "grid_util.h"
#include "util.h"
#include "staggered_grid.h"
#include "constants.h"
#include "util.h"

#include <igl/grid.h>
#include <igl/bounding_box.h>

#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <stdio.h>
#include <time.h>
#include <iostream>

typedef Eigen::Triplet<double> T;

// Predefined colors
const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
const Eigen::RowVector3d yellow(1.0, 0.9, 0.2);
const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
const Eigen::RowVector3d green(0.2, 0.6, 0.3);
const Eigen::RowVector3d black(0.0, 0.0, 0.0);
const Eigen::RowVector3d white(1.0, 1.0, 1.0);
const Eigen::RowVector3d red(0.8, 0.2, 0.2);

// Helper variables
StaggeredGrid staggeredGrid;

// Viewer data ids
int smokeId;
int boxId;
int ballId;
int bunnyId;

// == Simulation State ==
Eigen::MatrixXd q;

Eigen::RowVector3d currBallCenter;
Eigen::MatrixXd ballV;
Eigen::MatrixXi ballF;

Eigen::RowVector3d currBunnyCenter;
Eigen::MatrixXd bunnyV;
Eigen::MatrixXi bunnyF;

Eigen::MatrixXd particleTemplateV;
Eigen::MatrixXi particleTemplateF;
Eigen::MatrixXd particlesV;
Eigen::MatrixXi particlesF;
// ====

Eigen::RowVector3d bunnyBoxHalfLengths;

// helps convert points to particles
Eigen::SparseMatrix<double> selectionMatrix;
Eigen::MatrixXd templateParticlesV;

// Update location and velocity of smoke particles
inline void simulate()
{
	staggeredGrid.updateExternalObjects(currBallCenter, &bunnyV, &bunnyF);

	// 1. update velocities
	staggeredGrid.advectVelocities();
	staggeredGrid.applyExternalForces();
	staggeredGrid.applyPressureProjections();

	// 2. advect temperature and density
	staggeredGrid.updateTemperatureAndDensity();

	// 3. advect particles
	staggeredGrid.advectPosition(q);
}

inline void simulateBall()
{
	Eigen::RowVector3d ballMaxCorner = currBallCenter + Eigen::RowVector3d::Constant(ballRadius) + dt * ballVelocity;
	Eigen::RowVector3d ballMinCorner = currBallCenter - Eigen::RowVector3d::Constant(ballRadius) + dt * ballVelocity;
	if (isInBox(SMOKE_BOX, ballMinCorner) && isInBox(SMOKE_BOX, ballMaxCorner))
	{
		for (int i = 0; i < ballV.rows(); i++)
		{
			ballV.row(i) += dt * ballVelocity;
		}
		currBallCenter += dt * ballVelocity;
	}
}

inline void simulateBunny()
{
	Eigen::RowVector3d bunnyMaxCorner = currBunnyCenter + bunnyBoxHalfLengths + dt * bunnyVelocity;
	Eigen::RowVector3d bunnyMinCorner = currBunnyCenter - bunnyBoxHalfLengths + dt * bunnyVelocity;
	if (isInBox(SMOKE_BOX, bunnyMaxCorner) && isInBox(SMOKE_BOX, bunnyMinCorner))
	{
		for (int i = 0; i < bunnyV.rows(); i++)
		{
			bunnyV.row(i) += dt * bunnyVelocity;
		}
		currBunnyCenter += dt * bunnyVelocity;
	}
}

inline void createSmokeBox(Eigen::MatrixXd& boxV, Eigen::MatrixXi& boxF, Eigen::MatrixXd& q, const Eigen::AlignedBox3d& boundary)
{
	// Create box
	igl::read_triangle_mesh("../data/box.obj", boxV, boxF);

	transformVertices(boxV, boundary);

	// Create smoke particles inside box
	q.resize(PARTICLE_COUNT, 3);
	Eigen::Vector3d bottomLeftFloor = SMOKE_BOUNDS.corner(Eigen::AlignedBox3d::BottomLeftFloor);
	Eigen::Vector3d topRightCeil = SMOKE_BOUNDS.corner(Eigen::AlignedBox3d::TopRightCeil);
	for (int i = 0; i < PARTICLE_COUNT; i++)
	{
		q(i, 0) = getRand(bottomLeftFloor(0), topRightCeil(0));
		q(i, 1) = getRand(bottomLeftFloor(1), topRightCeil(1));
		q(i, 2) = getRand(bottomLeftFloor(2), topRightCeil(2));
	}
}

inline void initializeBall()
{
   	igl::read_triangle_mesh("../data/sphere.obj", ballV, ballF);

   	// Display sphere smaller than actual size to account for the particle's large size
   	Eigen::AlignedBox3d sphereBoundaries;
   	sphereBoundaries.extend(initialBallPosition - Eigen::Vector3d::Constant(ballRadius - 0.5));
   	sphereBoundaries.extend(initialBallPosition + Eigen::Vector3d::Constant(ballRadius - 0.5));
   	transformVertices(ballV, sphereBoundaries);
}

inline void initializeBunny()
{
	igl::read_triangle_mesh("../data/coarser_bunny.obj", bunnyV, bunnyF);

	// Find dimensions of bunny
	Eigen::MatrixXd boundingV;
	Eigen::MatrixXi boundingF;
	igl::bounding_box(bunnyV, boundingV, boundingF);
	Eigen::AlignedBox3d defaultBox;
	createAlignedBox(boundingV, defaultBox);

	// Scale bunny to set size
	double scale = bunnyHalfLength / (defaultBox.sizes()(0) / 2.0);
	bunnyBoxHalfLengths = Eigen::RowVector3d(bunnyHalfLength, defaultBox.sizes()(1) * scale / 2, defaultBox.sizes()(2) * scale / 2);
	Eigen::Vector3d minCorner = initialBunnyPosition - bunnyBoxHalfLengths.transpose();
	Eigen::Vector3d maxCorner = initialBunnyPosition + bunnyBoxHalfLengths.transpose();
	Eigen::AlignedBox3d newBox(minCorner, maxCorner);

	transformVertices(bunnyV, newBox);
}

inline void initializeParticles()
{
	// Create template particle
	igl::readOBJ("../data/particle.obj", particleTemplateV, particleTemplateF);
	transformVertices(particleTemplateV, PARTICLE_SIZE);

	int stepV = particleTemplateV.rows();
	int stepF = particleTemplateF.rows();

	particlesV.resize(PARTICLE_COUNT * stepV, 3);
	particlesV = particleTemplateV.replicate(PARTICLE_COUNT, 1);

	Eigen::MatrixXi ones;
	ones.resize(particleTemplateF.rows(), particleTemplateF.cols());
	ones.setOnes();
	particlesF.resize(PARTICLE_COUNT * stepF, 3);
	for (int i = 0; i < PARTICLE_COUNT; i++)
	{
		particlesF.block(i * stepF, 0, stepF, 3) = particleTemplateF + ones * (i * stepV);
	}

	// initialize selection matrix to convert points to particles
	std::vector<T> sparseEntries;
	for (int i = 0; i < PARTICLE_COUNT; i++)
	{
		for (int j = 0; j < stepV; j++)
		{
			sparseEntries.push_back(T(i * stepV + j, i, 1));
		}
	}
	selectionMatrix.resize(PARTICLE_COUNT * stepV, PARTICLE_COUNT);
	selectionMatrix.setZero();
	selectionMatrix.setFromTriplets(sparseEntries.begin(), sparseEntries.end());

	templateParticlesV = particleTemplateV.replicate(PARTICLE_COUNT, 1);
}

inline void updateParticleMeshes()
{
	int stepV = particleTemplateV.rows();
	//for (int i = 0; i < q.rows(); i++)
	//{
	//	// translate template particle to appropriate location
	//	particlesV.block(i * stepV, 0, stepV, 3) 
	//		= particleTemplateV + q.row(i).replicate(stepV, 1);
	//}
	particlesV = templateParticlesV + selectionMatrix * q;
}

// Must be called first
inline void simulation_setup(int argc, char** argv)
{
	// Add sphere
	if (ball)
	{
		currBallCenter = initialBallPosition.transpose();
		initializeBall();
		ballId = Visualize::addObjectToScene(ballV, ballF, orange);
	}

	// Add bunny
	if (bunny)
	{
		currBunnyCenter = initialBunnyPosition.transpose();
		initializeBunny();
		bunnyId = Visualize::addObjectToScene(bunnyV, bunnyF, orange);
	}

	// Add smoke and box
	Eigen::MatrixXd boxV;
	Eigen::MatrixXi boxF;
	createSmokeBox(boxV, boxF, q, SMOKE_BOX);

	if (useParticles)
	{
		// Create all particles
		initializeParticles();
		updateParticleMeshes();
		smokeId = Visualize::addObjectToScene(particlesV, particlesF, white);
		Visualize::viewer().data(smokeId).set_face_based(false);
		Visualize::viewer().core().toggle(Visualize::viewer().data(smokeId).show_lines);
	}
	else
	{
		smokeId = Visualize::addPointsToScene(q, white);
	}

	// Create box last to set camera focus on the box
	boxId = Visualize::addObjectToScene(boxV, boxF, orange);
	Visualize::setInvisible(boxId, true);

	staggeredGrid = StaggeredGrid(q, SMOKE_BOX, GRID_DIM, SMOKE_BOX.sizes()(0) / (GRID_DIM(0) - 1.0));

	////// TODO: DELETE. Testing if initialization of staggered grid points is correct
	//Eigen::MatrixXd u, v, w, p;
	//staggeredGrid.getGridPoints(u, v, w, p);
	//Visualize::addPointsToScene(u, yellow);
	//Visualize::addPointsToScene(v, orange);
	//Visualize::addPointsToScene(w, green);
	//Visualize::addPointsToScene(p, red);
}

inline void draw()
{
	if (useParticles)
	{
		updateParticleMeshes();
		Visualize::updateObject(smokeId, particlesV);
	}
	else
	{
		Visualize::updatePoints(smokeId, q, white);
	}

	if (ball)
	{
		Visualize::updateObject(ballId, ballV);
	}
	
	if (bunny)
	{
		Visualize::updateObject(bunnyId, bunnyV);
	}
}

#endif
