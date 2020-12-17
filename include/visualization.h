#ifndef  VISUALIZATION_H
#define  VISUALIZATION_H

#include <Eigen/Dense>

#include <imgui/imgui.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

/*
* Wrapper around libigl's viewer and menu objects
*/
namespace Visualize {
	igl::opengl::glfw::Viewer& viewer();
	void setup(const Eigen::MatrixXd& q, const Eigen::MatrixXd& qdot);
	int addObjectToScene(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::RowVector3d& color);
	void setInvisible(int dataId);
	void addPointsToScene(int dataId, const Eigen::MatrixXd& points, const Eigen::RowVector3d& color);
}

#endif