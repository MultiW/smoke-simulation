#include <iostream>
#include <thread>

#include "visualization.h"
#include "simulation.h"

// Simulation state
// TODO: what format is required of q?
Eigen::MatrixXd q;
Eigen::MatrixXd qdot;

//simulation time and time step
double t = 0; //simulation time 
double dt = 0.001; //time step

bool simulation_callback()
{
	while (true)
	{
		simulate(q, qdot, dt, t);
		t += dt;
	}
	return false;
}

bool draw_callback(igl::opengl::glfw::Viewer &viewer) 
{
    draw(q, qdot, t);

    return false;
}

int main(int argc, char* argv[])
{
    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    //setup libigl viewer and activate 
    Visualize::setup(q, qdot);
    Visualize::viewer().callback_post_draw = &draw_callback;
    Visualize::viewer().launch();

    return 1; 
}
