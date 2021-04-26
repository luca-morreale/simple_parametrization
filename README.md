# Simple Parametrization

This repo containes implementations related to parametrization algorithms.

## Algorthims
- SLIM (with and without fixed boundary)
- Tutte

## Utils
- Mesh cut
- Seam path definition

## Requirements
All the algorithms rely on [CGAL](https://www.cgal.org/) or [libigl](https://libigl.github.io/).
Other dependencies:
- boost
- Eigen3
- gmp
- cmake

To install libigl is it enough to link it to inside the `include` folder:
```sh
mkdir include
cd include
ln -s /your/path/to/libigl/include/igl igl
ln -s /your/path/to/libigl/external/glad/include/glad glad
```

## Compile
```sh
mkdir build
cd build
cmake .. && make -j4
```

## SLIM
SLIM parametrization can be compute with:
```sh
./build/slim your_obj_file.obj
```
this will generate a new file `your_obj_file_slim.obj` containing the parametrize mesh to the unit disk.
Note this algorithm requires a surface homeomorphic to a disk. If has a higher genus, then it is possible to cut it.

## Free SLIM
SLIM parametrization without a fixed boundary can be compute with:
```sh
./build/slim your_obj_file.obj
```
this will generate a new file `your_obj_file_freeslim.obj` containing the parametrize mesh.
Note this algorithm requires a surface homeomorphic to a disk. If has a higher genus, then it is possible to cut it.

## Tutte
SLIM parametrization can be compute with:
```sh
./build/tutte your_obj_file.obj
```
this will generate a new file `your_obj_file_tutte.obj` containing the parametrize mesh to a disk.
Note this algorithm requires a surface homeomorphic to a disk. If has a higher genus, then it is possible to cut it.


## Seam path definition
A seam path is computed based on an ordered set of points it cut through.
These points can be defined in pairs. For example:
```
4199 3041
3041 410
```
the seam will list all vertex starting from vertex #4199 to vertex #3041. Then from vertex #3041 to vertex #410.
The seam path is defined based on Dijkstra algorithm.

Running:
```sh
./build/dijkstra_seam your_obj_file.off your_vertex_list.txt
```
produces a txt file `your_obj_file.off.selection.txt` containing the seam path.

## Seam cut
Once the seam path is defined, it is possible to cut the mesh:
```sh
./build/cut your_obj_file.off your_obj_file.off.selection.txt
```
this produces a new obj file `your_obj_file_cut.obj`.

