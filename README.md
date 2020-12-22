# Marker-and-Cell (MAC) Based 2D Fluid Simulation with `libigl`

## Summary

The final project for the course Physics Animations at University of Toronto. Currently only a 2D grid based simulation is available. ([final report](https://docs.google.com/document/d/1vo__Ta67T-i1MSSQ-jYLqVYwf1TI6slXclJ3VWIcJeo/edit?usp=sharing))

### Prerequisite installation

On all platforms, we will assume you have installed cmake and a modern c++
compiler on Mac OS X[¹](#¹macusers), Linux[²](#²linuxusers), or
Windows[³](#³windowsusers).

We also assume that you have cloned this repository using the `--recursive`
flag (if not then issue `git submodule update --init --recursive`).

**Note:** We only officially support these assignments on Ubuntu Linux 18.04 (the OS the teaching labs are running) and OSX 10.13 (the OS I use on my personal laptop). While they *should* work on other operating systems, we make no guarantees.

## Usage

1. clone the repository using `--recursive` flag to include `libigl` when cloning (if not then issue `git submodule update --init --recursive`).
2. run `mkdir build` to create the build directory, then run `cd build` to go into the directory
3. run `cmake .. -DCMAKE_BUILD_TYPE=Release` or `cmake .. -DCMAKE_BUILD_TYPE=Debug` for debug mode inside the `build` directory
4. run `make` inside the `build` directory
5. run `./fluid-sim` for the default simulation with a 50x50 2d grid

## Dependencies

Libraries used: Eigen, Libigl.

## References

- [SIGGRAPH 2007 Course Notes](https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf)
- [CSC417 Physics Based Animation](https://github.com/dilevin/CSC417-physics-based-animation)
- [Fluid Simulation Slides](https://github.com/dilevin/CSC417-physics-based-animation/blob/master/lectures/10-fluid-simulation-final.pdf)
