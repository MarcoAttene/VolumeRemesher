# Volume remeshing

This code implements a variation of the volume meshing algorithm described in "**Convex Polyhedral Meshing for Robust Solid Modeling**" by Lorenzo Diazzi and Marco Attene (ACM Trans Graphics Vol 40, N. 6, Procs of SIGGRAPH Asia 2021). You may download a copy of the paper preprint here: http://arxiv.org/abs/2109.14434
The original code is available here: https://github.com/MarcoAttene/VolumeMesher

In this version, a new function is made available to embed a set of triangles within an existing tetrahedral mesh.
See embed.h to start.

## Usage
Clone this repository, including submodules, with:
```
git clone --recursive https://github.com/MarcoAttene/VolumeRemesher
```
Once done, you may build the executable as follows:
```
mkdir build
cd build
cmake ..
```
This will produce an appropriate building configuration for your system.
On Windows MSVC, this will produce a mesh_generator.sln file.
On Linux/MacOS, this will produce a Makefile. 
Use it as usual to compile mesh_generator.

When compiled, the code generates an executable called ``mesh_generator``.
Launch it with no command line parameters to have a list of supported options.

Example:

```
mesh_generator input_triangles.off U input_tetmesh.tet
```
creates a file called ``volume.msh`` representing the (subdivided) volume mesh with input triangles embedded.

The input tetrahedral mesh may be represented in either TET or MEDIT format (both ASCII).

All the other switches and operational modes available in the original https://github.com/MarcoAttene/VolumeMesher are still available here.


We tested our code on MacOS (GCC-10) and Windows (MSVC 2019).
It should work on Linux-GCC and MacOS-Clang too, but we have not tested it on these configurations.

|:warning: WARNING: Apparently, CLANG does not support a fully IEEE compliant floating point environment which is necessary to guarantee that indirect predicates work as expected. The only way we found to guarantee correctness on this compiler was to disable all optimizations. Please be aware of this fact should you notice a performance degradation in your experiments. |
| --- |
