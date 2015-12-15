# SPH Fluid Simulation

### Introduction

This code was written for the semester project of the Physically Based Simulation course in 2015 at ETHZ. The goal was to implement a basic SPH fluid simulator in C++11. Ideas for the implementation are taken from the following papers:

- [1] [Weakly compressible SPH for free surface flows](http://cg.informatik.uni-freiburg.de/publications/2007_SCA_SPH.pdf)
- [2] [Predictive-Corrective Incompressible SPH](https://graphics.ethz.ch/~sobarbar/papers/Sol09/Sol09.pdf)
- [3] [Boundary Handling handling and adaptive time-stepping for PCISPH](http://cg.informatik.uni-freiburg.de/publications/2010_VRIPHYS_boundaryHandling.pdf)
- [4] [Versatile Rigid-Fluid Coupling for Incompressible SPH](http://cg.informatik.uni-freiburg.de/publications/2012_SIGGRAPH_rigidFluidCoupling.pdf)
- [5] [Versatile Surface Tension and Adhesion for SPH Fluids](http://cg.informatik.uni-freiburg.de/publications/2013_SIGGRAPHASIA_surfaceTensionAdhesion.pdf)
- [6] [Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels](http://www.cc.gatech.edu/~turk/my_papers/sph_surfaces.pdf)

### Features

The following features have been implemented:

- Basic OpenGL visualization
- JSON based scene description
- Index-sorted uniform grid for neighbour search
- Wavefront OBJ support
- WCSPH [1] and PCISPH [2] solvers
- PCISPH with adaptive time-stepping [3]
- Boundaries using boundary particles [3], [4]
    - Create boundary particles for boxes, spheres and arbitrary meshes
- Surface tension forces [5]
- Fluid mesh generation using marching cubes
    - Isotropic kernel
    - Anisotropic kernel (only works partially yet) [6]
- Cache system to cache particles and meshes
- Render application using either OpenGL visualization or SmallLuxGPU4

### Results

TODO

### Third Party Code

Part of the framework was derived from the [Nori](https://github.com/wjakob/nori) codebase written by Wenzel Jakob. In addition, the following libraries have been used:

- [cxxopts](https://github.com/jarro2783/cxxopts)
- [Eigen](http://eigen.tuxfamily.org)
- [filesystem](https://github.com/wjakob/filesystem)
- [json11](https://github.com/dropbox/json11)
- [libexecstream](http://libexecstream.sourceforge.net)
- [nanogui](https://github.com/wjakob/nanogui)
- [pcg32](https://github.com/wjakob/pcg32)
- [stb](https://github.com/nothings/stb)
- [tbb](https://www.threadingbuildingblocks.org)
- [tinydir](https://github.com/cxong/tinydir)
- [tinyformat](https://github.com/c42f/tinyformat)

