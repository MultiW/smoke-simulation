#include "visualization.h"
#include "constants.h"
#include "grid_util.h"
#include "util.h"

#include <igl/read_triangle_mesh.h>

#include <stdio.h>
#include <iostream>

namespace Visualize
{
	igl::opengl::glfw::Viewer g_viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	bool firstDataUsed = false; // has first data been used

	Eigen::MatrixXd particlesV;
	Eigen::MatrixXi particlesF;
}

igl::opengl::glfw::Viewer& Visualize::viewer()
{
	return g_viewer;
}

void draw_main_viewer_menu()
{
	// Draw parent menu content
	Visualize::menu.draw_viewer_menu();
}

void Visualize::setup()
{
	// Set up menu
	g_viewer.plugins.push_back(&menu);
	menu.callback_draw_viewer_menu = draw_main_viewer_menu;
}

int Visualize::addPointsToScene(const Eigen::MatrixXd& points, const Eigen::RowVector3d& color)
{
	int dataId;
	if (firstDataUsed)
	{
		g_viewer.append_mesh();
	}
	firstDataUsed = true;
	dataId = g_viewer.data().id;

	g_viewer.data(dataId).set_points(points, color);
	return dataId;
}

int Visualize::addObjectToScene(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::RowVector3d& color)
{
	if (firstDataUsed)
	{
		g_viewer.append_mesh();
	}
	firstDataUsed = true;
	int dataId = g_viewer.data().id;

	g_viewer.data(dataId).set_mesh(V, F);
	g_viewer.data(dataId).set_colors(color);
	g_viewer.data(dataId).set_face_based(true);
	return dataId;

}

void Visualize::updatePoints(int dataId, const Eigen::MatrixXd& V, const Eigen::RowVector3d& color)
{
	g_viewer.data(dataId).set_points(V, color);
}

void Visualize::updateObject(int dataId, const Eigen::MatrixXd& V)
{
	g_viewer.data(dataId).set_vertices(V);
}

void Visualize::setInvisible(int dataId, bool status)
{
	g_viewer.data(dataId).show_faces = !status;
	// TODO: hide lines too
}

