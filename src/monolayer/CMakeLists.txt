
include_directories(.)

add_library(planar_geometry STATIC GraphConstructor.cpp PlanarGeometry.cpp Utils.cpp)

# Executable to get area and perimeter statistics of film
add_executable(extractFilmParams extractFilmParams.cpp)
target_link_libraries(extractFilmParams PUBLIC planar_geometry OpenMP::OpenMP_CXX ${GROMACS_LIBRARIES})

# Executable to get 2Drdf distribution from one trajectory with both kinds of particles
add_executable(extractPlanarDistribution extractPlanarDistribution.cpp)
target_link_libraries(extractPlanarDistribution PUBLIC planar_geometry OpenMP::OpenMP_CXX ${GROMACS_LIBRARIES})

# Test for the graph2d library
add_executable(GraphTest GraphTest.cpp)
target_link_libraries(GraphTest PUBLIC planar_geometry OpenMP::OpenMP_CXX)

# Test for the graph2d library with 400K data
add_executable(GraphTest400K GraphTest400K.cpp)
target_link_libraries(GraphTest400K PUBLIC planar_geometry OpenMP::OpenMP_CXX)
