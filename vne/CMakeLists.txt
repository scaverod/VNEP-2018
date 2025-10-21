cmake_minimum_required(VERSION 2.7)
project(vne)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(BOOST_DIR /usr/local/boost/)
set(Boost_INCLUDE_DIR /usr/local/boost)
set(Boost_LIBRARY_DIR /usr/local/boost/libs)
FIND_PACKAGE(Boost 1.54.0 REQUIRED)
FIND_PACKAGE(cplex REQUIRED)

# Debug output
message(STATUS "CPLEX_INCLUDE_DIR: ${CPLEX_INCLUDE_DIR}")
message(STATUS "CONCERT_INCLUDE_DIR: ${CONCERT_INCLUDE_DIR}")
message(STATUS "CPLEX_INCLUDE_DIRS: ${CPLEX_INCLUDE_DIRS}")
message(STATUS "CPLEX_LIBRARY: ${CPLEX_LIBRARY}")
message(STATUS "CPLEX_LIBRARIES: ${CPLEX_LIBRARIES}")

#find_package(Boost COMPONENTS program_options timer system chrono)
# Use CPLEX_INCLUDE_DIRS if available, otherwise fall back to individual directories
if(CPLEX_INCLUDE_DIRS)
  include_directories(include ${Boost_INCLUDE_DIR} ${CPLEX_INCLUDE_DIRS})
else()
  include_directories(include ${Boost_INCLUDE_DIR} ${CPLEX_INCLUDE_DIR} ${CONCERT_INCLUDE_DIR})
endif()
#link_directories(${Boost_LIBRARY_DIR} ${CPLEX_LIBRARIES})

# set(cplex_DIR /opt/cplex)


set(CMAKE_CXX_FLAGS "-Wno-unknown-pragmas -g -pedantic -lm -std=c++0x -O3 -fPIC -DIL_STD -fexceptions -m64")

set(SRC src/main.cpp src/backtrack.cpp  src/bb.cpp src/bbcgu.cpp src/bprice_rec.cpp src/bprice.cpp  src/cg.cpp  src/cgunsplittable.cpp  src/cplex.cpp  src/fix.cpp src/ac.cpp  src/util.cpp  src/vnedefs.cpp)
ADD_EXECUTABLE(vne.out ${SRC})
TARGET_LINK_LIBRARIES(vne.out ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} dl)
#message(${BOOST_library_dirs})

# tests
#add_executable(test.out test/main.cpp ${SRC})
#TARGET_LINK_LIBRARIES(test.out ${CPLEX_LIBRARIES} ${Boost_LIBRARIES})

#enable_testing()

#add_test(backtrackXcplex test.out ../data/simple/sub1test.gtt ../data/simple/vir3test.gtt)
