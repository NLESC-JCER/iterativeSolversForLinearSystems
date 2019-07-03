# Install script for directory: /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/AdolcForward"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/AlignedVector3"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/ArpackSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/AutoDiff"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/BVH"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/EulerAngles"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/FFT"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/IterativeSolvers"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/KroneckerProduct"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/LevenbergMarquardt"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/MatrixFunctions"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/MoreVectorization"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/MPRealSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/NonLinearOptimization"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/NumericalDiff"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/OpenGLSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/Polynomials"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/Skyline"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/SparseExtra"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/SpecialFunctions"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

