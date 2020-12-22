# MAC based fluid simulation based on the SIGGRAPH 2007 Course Notes from Robert Bridson
[Link to notes](https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf)

## Summary
The final project for the course Physics Animations at University of Toronto
Currently only a 2D grid based simulation is available.


## Implementation
Libraries used: Eigen, Libigl.



## Compilation and usage
1. clone the repository using `--recursive` flag to include `libigl` when cloning.
2. run `mkdir build` to create the build directory, then run `cd build` to go into the directory
3. run `cmake .. -DCMAKE_BUILD_TYPE=Release` or `cmake .. -DCMAKE_BUILD_TYPE=Debug` for debug mode inside the `build` directory
4. run `make` inside the `build` directory
5. run `./fluid-sim` for the default simulation with a 50x50 2d grid


