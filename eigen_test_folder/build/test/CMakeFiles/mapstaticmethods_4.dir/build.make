# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build

# Include any dependencies generated for this target.
include test/CMakeFiles/mapstaticmethods_4.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/mapstaticmethods_4.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/mapstaticmethods_4.dir/flags.make

test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o: test/CMakeFiles/mapstaticmethods_4.dir/flags.make
test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o: ../test/mapstaticmethods.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o"
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o -c /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/test/mapstaticmethods.cpp

test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.i"
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/test/mapstaticmethods.cpp > CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.i

test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.s"
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/test/mapstaticmethods.cpp -o CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.s

test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o.requires:

.PHONY : test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o.requires

test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o.provides: test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/mapstaticmethods_4.dir/build.make test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o.provides.build
.PHONY : test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o.provides

test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o.provides.build: test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o


# Object files for target mapstaticmethods_4
mapstaticmethods_4_OBJECTS = \
"CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o"

# External object files for target mapstaticmethods_4
mapstaticmethods_4_EXTERNAL_OBJECTS =

test/mapstaticmethods_4: test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o
test/mapstaticmethods_4: test/CMakeFiles/mapstaticmethods_4.dir/build.make
test/mapstaticmethods_4: test/CMakeFiles/mapstaticmethods_4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mapstaticmethods_4"
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mapstaticmethods_4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/mapstaticmethods_4.dir/build: test/mapstaticmethods_4

.PHONY : test/CMakeFiles/mapstaticmethods_4.dir/build

test/CMakeFiles/mapstaticmethods_4.dir/requires: test/CMakeFiles/mapstaticmethods_4.dir/mapstaticmethods.cpp.o.requires

.PHONY : test/CMakeFiles/mapstaticmethods_4.dir/requires

test/CMakeFiles/mapstaticmethods_4.dir/clean:
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/test && $(CMAKE_COMMAND) -P CMakeFiles/mapstaticmethods_4.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/mapstaticmethods_4.dir/clean

test/CMakeFiles/mapstaticmethods_4.dir/depend:
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/test /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/test /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/test/CMakeFiles/mapstaticmethods_4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/mapstaticmethods_4.dir/depend

