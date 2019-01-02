# Toy-Fluid-Simulation
A simple fluid simulator using physics-based method guided by these two tutorials: http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf and https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf.

- Dependencies:
	- C++11

	- OpenGL core profile version 3.3 (with [glfw](https://www.glfw.org/download.html), [glad](https://glad.dav1d.de/) downloaded)
	
	- [Eigen 3.3](http://eigen.tuxfamily.org/)
	
	-  gcc compiler (recommended)
	
- Compile instructions: 
In order to compile and execute the project correctly, you'll want to have the requirements above filled. Pay attention to the [glad](https://glad.dav1d.de/) item, you **SHOULD** download the right `glad.c` file from the [website]((https://glad.dav1d.de/)) according to your specific environment and modify `PATH_TO_FILE_glad.c` in `CMakeList.txt` with the correct path to it. 
The remaining part is quite simple. Suppose you are now in the directory containing the `CMakeList.txt` file and all those source codes, just build the project with the following commands:
	- `mkdir build`
	- `cd build`
	- `cmake ..`
	- `make` or `make -jX` (where X is the number of CPU cores you want to use)
	
	Hopefully the executable file named `ToyFluidSimProj` appears in the build directory. Click it and it will start simulating until `Esc` key is pressed.
 
- Notes:
	- This project offers the flexibility to set some user-defined parameters: number of grid cells in x-direction (`Nx`) and in y-direction (`Ny`), position of the inlets, fluid velocity at the inlet, and viscosity of the fluid(`visc`). `Nx` and `Ny` are global const variables at the top of the `main.cpp` file, while parameters concerning the inlets is set at the beginning of the main function (see `smlt.setInlet()`), and `visc` is set through `smlt.setVisc()`.
	
	- Some examples of simulation can be found in the `examples` directory. Their corresponding settings are listed below:
		- fluid1: Nx=200, Ny=100, smlt.setInlet(3, {33, 33}, {3, 3}, 3, {33, 66}, {3, -3}), no viscosity
		
		- fluid2: Nx=200, Ny=100, smlt.setInlet(2, {33, 66}, {3, 0}, 2, {166, 66}, {-3, 0}), no viscosity
		
		- fluid3: Nx=200, Ny=100, smlt.setInlet(2, {33, 66}, {3, 0}, 2, {166, 66}, {-3, 0}), visc=0.07
		
		- fluid4: Nx=200, Ny=100, smlt.setInlet(2, {33, 66}, {3, 0}, 2, {166, 66}, {-3, 0}), visc=0.3








