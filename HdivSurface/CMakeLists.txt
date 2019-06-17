# This project works with H(div) on manifolds and MHM formulation.

set(common_SRC
main.cpp
Tools.h
Tools.cpp
ProblemConfig.h
)

add_executable(Surface ${common_SRC} main.cpp)

target_link_libraries(Surface pz)



configure_file(../Meshes/Quad.msh Quad.msh COPYONLY)
configure_file(../Meshes/Cube.msh Cube.msh COPYONLY)
configure_file(../Meshes/BasicMesh.msh BasicMesh.msh COPYONLY)
configure_file(../Meshes/LMesh.msh LMesh.msh COPYONLY)
configure_file(../Meshes/esfera.msh esfera.msh COPYONLY)
configure_file(../Meshes/esfera2.msh esfera.msh COPYONLY)
configure_file(../Meshes/Cube.msh esfera.msh COPYONLY)
