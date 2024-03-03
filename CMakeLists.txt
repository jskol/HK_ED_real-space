  cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)
	
  project(myproject LANGUAGES CXX)
	

  if(NOT CMAKE_BUILD_TYPE)
  	set(CMAKE_BUILD_TYPE Release)
  endif() 
  
  set(CMAKE_CXX_FLAGS "-Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")

  message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

  set(CMAKE_CXX_STANDARD 17)
  #set(libcommute_DIR ${libcommute_ROOT}/C:/Users/avoga/Documents/libcommute.src/cmake)
  set(libcommute_DIR /home/krystian/Dokumenty/libcommute.build/home/Dokumenty/lib/cmake/libcommute)

  #set (CMAKE_PREFIX_PATH ="C:/Users/avoga/Documents/libcommute.src/cmake")
  #set (CMAKE_PREFIX_PATH "C:/Users/avoga/Documents/libcommute.build")
  find_package(OpenMP)
  
  # Import libcommute target
  find_package(libcommute 0.7.2 CONFIG REQUIRED)
  
  # Import lambda lanczos
  add_library(lambda_lanczos INTERFACE)
  target_include_directories(lambda_lanczos INTERFACE /home/krystian/Dokumenty/lambda-lanczos/include/lambda_lanczos)

  # Import eigen
  #set(Eigen_DIR $ENV{HOME}/Tools/eigen/cmake)
  set (CMAKE_PREFIX_PATH "C:/Users/avoga/Documents/eigen3.build")
  find_package (Eigen3 3.3 REQUIRED NO_MODULE)  

  # Build an executable called 'H-K' 
  set( CODES HK-Edge_states)
  foreach(e ${CODES})
	  add_executable( ${e} ${e}.cpp)
  	  if(OpenMP_CXX_FOUND)
        message( "\n Compiling with OpenMP")
		    target_link_libraries(${e} PRIVATE libcommute::libcommute OpenMP::OpenMP_CXX Eigen3::Eigen lambda_lanczos)
	  else()
	    message("\n Compiling without OpenMP")
      target_link_libraries(${e} PRIVATE libcommute::libcommute Eigen3::Eigen lambda_lanczos)
  	  
  	  endif()
  endforeach()
	  
  
