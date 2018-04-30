find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Eigen3 REQUIRED)
add_library(eigen INTERFACE IMPORTED)

set(SOURCE_FILES main.cpp
		 simulator.cpp
		 smoke_sim.cpp
		 fluid_sim.cpp
		 mac_grid.cpp
		 util/vec.cpp
		 grid_data.cpp
		 camera.cpp
		 util/fps.cpp
		 util/constants.cpp
		 util/basic_math.cpp)
add_SMOKE_executable(SMOKE ${SOURCE_FILES})
include_directories( ${OPENGL_INCLUDE_DIR}  ${GLUT_INCLUDE_DIRS} )
target_link_libraries(SMOKE ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )
target_include_directories(SMOKE SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(SMOKE eigen)
target_link_libraries(SMOKE partio)