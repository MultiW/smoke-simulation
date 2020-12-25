#include <iostream>
#include <thread>

#include "visualization.h"
#include "simulation.h"
#include "constants.h"

// Simulation state
// TODO: what format is required of q?
Eigen::MatrixXd q;

//simulation time and time step
double t = 0; //simulation time 

bool simulation_callback()
{
	while (true)
	{
		simulate(q, t);
		t += dt;
	}
	return false;
}

bool draw_callback(igl::opengl::glfw::Viewer& viewer)
{
	draw(q, t);

	return false;
}

int main(int argc, char* argv[])
{
	simulation_setup(argc, argv, q);

	//run simulation in seperate thread to avoid slowing down the UI
	std::thread simulation_thread(simulation_callback);
	simulation_thread.detach();

	//setup libigl viewer and activate 
	Visualize::setup(q);
	Visualize::viewer().callback_post_draw = &draw_callback;
	Visualize::viewer().launch();

	return 1;
}
