CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( heat )
find_package(MPI)
set(CMAKE_CXX_COMPILER "mpic++")

include_directories(${MPI_INCLUDE_PATH})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -pedantic -Wall -Wextra -Wshadow -std=c++11")

option(ARCHITECTURE_OPTIMIZATION [ON])
option(ENABLE WERROR [ON])

SET (HEADERS Timer.h misc.h grid.h)

SET( SRCS
   cg.cpp grid.cpp solver.cpp
   # more *.cpp files follow here
   )
   
ADD_EXECUTABLE( heat ${HEADERS} ${SRCS})