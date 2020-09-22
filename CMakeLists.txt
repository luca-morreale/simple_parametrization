project(OptimalViewEstimator)
cmake_minimum_required(VERSION 2.8.3)

# ----------------------- LIBRARIES ---------------------------

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/)
set(CMAKE_BUILD_TYPE "Release")

# set(Boost_USE_STATIC_LIBS ON)
# set(Boost_USE_MULTITHREADED ON)
set(BOOST_ROOT "/usr/include/boost")

find_library(GMP_LIBRARY gmp /usr/lib)
find_package(Boost 1.65 REQUIRED COMPONENTS system thread)
find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# ----------------------- INCLUDES ----------------------------

include_directories(${EIGEN3_INCLUDE_DIR}
                      ${GMPXX_INCLUDE_DIR}
                      ${CGAL_INCLUDE_DIRS}
                      )
include_directories(${Boost_INCLUDE_DIR})
include_directories(include/)

# ----------------------- GCC FLAGS ----------------------------

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

# ------------------------ BUILD -------------------------------


add_executable(parametrize src/main_parametrize.cpp)
target_link_libraries(parametrize ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_LIBRARY} ${Boost_LIBRARIES} )


add_executable(cut src/main_cut.cpp)
target_link_libraries(cut ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_LIBRARY} ${Boost_LIBRARIES} )


