# Fluid Simulation: Smoke and Collision

## Demo

TODO: pictures

[Demonstration and methodology](https://youtu.be/lMkeVszKhB4)

## Code Navigation


## Simulation Options


## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git

C++11 or higher should be used.

## Compile
_Note: this project was only tested on Windows 10_

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `smoke-simulation.exe` binary.

More help with compilation can be found [here](http://libigl.github.io/libigl/tutorial/).

## Run

From within the `build` directory just issue:

    ./smoke-simulation.exe

A glfw app should launch displaying the simulation.

## Sources

TODO