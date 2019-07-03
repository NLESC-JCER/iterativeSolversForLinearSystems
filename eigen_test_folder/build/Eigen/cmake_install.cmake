# Install script for directory: /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Cholesky"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/CholmodSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Core"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Dense"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Eigen"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Eigenvalues"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Geometry"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Householder"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/IterativeLinearSolvers"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Jacobi"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/KLUSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/LU"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/MetisSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/OrderingMethods"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/PaStiXSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/PardisoSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/QR"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/QtAlignedMalloc"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/SPQRSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/SVD"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/Sparse"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/SparseCholesky"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/SparseCore"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/SparseLU"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/SparseQR"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/StdDeque"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/StdList"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/StdVector"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/SuperLUSupport"
    "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

