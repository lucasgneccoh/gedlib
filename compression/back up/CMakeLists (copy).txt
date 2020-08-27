cmake_minimum_required(VERSION 2.6)

project(COMPRESSION)
set(VERSION 1.0)

#set(EXECUTABLE_OUTPUT_PATH /bin})

# Define variables.

# Set up the compiler.

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic -Wredundant-decls -Wdisabled-optimization -Wctor-dtor-privacy -Wnon-virtual-dtor -Werror=old-style-cast -Wsign-promo -Werror=return-type -Wno-unused-parameter")
if(APPLE)
  set(CMAKE_CXX_COMPILER clang++ CACHE STRING "C++ compiler" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()


set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_MACOSX_RPATH ON)

set(EXEC_NAME arborescence_based_compression)

#set(GEDLIB /home/lucas/Documents/stage/gedlib)
set(GEDLIB ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(LIBRARY_OUTPUT_PATH ${GEDLIB}/lib)
set(LSAPE_HOME ${GEDLIB}/ext/lsape.5)
set(NOMAD_HOME ${GEDLIB}/ext/nomad.3.8.1)
set(LIBSVM_HOME ${GEDLIB}/ext/libsvm.3.22)
set(EIGEN_HOME ${GEDLIB}/ext/eigen.3.3.4)
set(FANN_HOME ${GEDLIB}/ext/fann.2.2.0)
set(BOOST_HOME ${GEDLIB}/ext/boost.1.69.0)

# Add include and link directories.

include_directories(SYSTEM ${GEDLIB} ${BOOST_HOME} ${NOMAD_HOME}/ext/sgtelib/src ${NOMAD_HOME}/src ${LSAPE_HOME}/include ${EIGEN_HOME}/Eigen ${LIBSVM_HOME} ${FANN_HOME}/include ${GUROBI_HOME}/include ${OMP_HOME}/include)


link_directories(${NOMAD_HOME}/lib ${FANN_HOME}/lib ${LIBSVM_HOME} ${GUROBI_HOME}/lib ${OMP_HOME}/lib)
#SHARED LIBRARY
link_directories(${GEDLIB}/lib)

#add_library(gxlgedlib SHARED ${GEDLIB}/src/env/ged_env.gxl.cpp)
#set_target_properties(gxlgedlib PROPERTIES SUFFIX ".so")
#target_link_libraries(gxlgedlib nomad doublefann svm)

#Génération de la liste des fichiers sources

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Add targets.

#add_custom_target(compression DEPENDS median_tests ibd_tests aids_edit_iso_test time_limit_tests)
#add_executable(compression src/compression_test.cpp)
#target_link_libraries(compression gxlgedlib)

add_subdirectory(src)




