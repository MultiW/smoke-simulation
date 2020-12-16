#ifndef  VISUALIZATION_H
#define  VISUALIZATION_H

#include <Eigen/Dense>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

/*
* Wrapper around libigl's viewer and menu objects
*/
namespace Visualize {
	igl::opengl::glfw::Viewer& viewer();
	void setup(const Eigen::VectorXd& q, const Eigen::VectorXd& qdot);
	int addObjectToScene(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::RowVector3d& color);
	void setInvisible(int dataId);
	void addPointsToScene(int dataId, const Eigen::MatrixXd& points, const Eigen::RowVector3d& color);
}

#endif