#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <Eigen/Dense>

#include <imgui/imgui.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

/*
* Wrapper around libigl's viewer and menu objects
*/
namespace Visualize {
	void setup();

	int addObjectToScene(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::RowVector3d& color);
	int addPointsToScene(const Eigen::MatrixXd& points, const Eigen::RowVector3d& color);
	int addParticlesToScene(const Eigen::MatrixXd& points, const Eigen::RowVector3d& color);

	void updatePoints(int dataId, const Eigen::MatrixXd& V, const Eigen::RowVector3d& color);
	void updateObject(int dataId, const Eigen::MatrixXd& V);
	void updateParticle(int dataId, const Eigen::MatrixXd& V);

	void setInvisible(int dataId, bool status);

	igl::opengl::glfw::Viewer& viewer();
}

#endif