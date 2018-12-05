cmake_minimum_required(VERSION 3.0)
project(kalmanfilter)

find_package(
 catkin_simple REQUIRED
 Armadillo REQUIRED)

catkin_simple()

cs_add_executable(${PROJECT_NAME} src/kalman_filter.cpp)
target_link_libraries(kalmanfilter ${ARMADILLO_LIBRARIES})

cs_install()

cs_export()
