#include <iostream>
#include <thread>

#include "visualization.h"
#include "simulation.h"
#include "constants.h"

//simulation time and time step
double t = 0; //simulation time 

bool simulation_callback()
{
	while (true)
	{
		// Move particles
		simulate();

		if (ball)
		{
			// Move sphere
			simulateBall();
		}

		t += dt;
	}
	return false;
}

bool draw_callback(igl::opengl::glfw::Viewer& viewer)
{
	draw();

	return false;
}

int main(int argc, char* argv[])
{
	simulation_setup(argc, argv);

	//run simulation in seperate thread to avoid slowing down the UI
	std::thread simulation_thread(simulation_callback);
	simulation_thread.detach();

	//setup libigl viewer and activate 
	Visualize::setup();
	Visualize::viewer().callback_post_draw = &draw_callback;
	Visualize::viewer().launch();

	return 1;
}
