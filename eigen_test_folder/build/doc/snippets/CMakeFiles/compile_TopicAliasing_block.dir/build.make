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
include doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/depend.make

# Include the progress variables for this target.
include doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/progress.make

# Include the compile flags for this target's objects.
include doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/flags.make

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o: doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/flags.make
doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o: doc/snippets/compile_TopicAliasing_block.cpp
doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o: ../doc/snippets/TopicAliasing_block.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o"
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o -c /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets/compile_TopicAliasing_block.cpp

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.i"
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets/compile_TopicAliasing_block.cpp > CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.i

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.s"
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets/compile_TopicAliasing_block.cpp -o CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.s

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o.requires:

.PHONY : doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o.requires

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o.provides: doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o.requires
	$(MAKE) -f doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/build.make doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o.provides.build
.PHONY : doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o.provides

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o.provides.build: doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o


# Object files for target compile_TopicAliasing_block
compile_TopicAliasing_block_OBJECTS = \
"CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o"

# External object files for target compile_TopicAliasing_block
compile_TopicAliasing_block_EXTERNAL_OBJECTS =

doc/snippets/compile_TopicAliasing_block: doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o
doc/snippets/compile_TopicAliasing_block: doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/build.make
doc/snippets/compile_TopicAliasing_block: doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable compile_TopicAliasing_block"
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_TopicAliasing_block.dir/link.txt --verbose=$(VERBOSE)
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets && ./compile_TopicAliasing_block >/home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets/TopicAliasing_block.out

# Rule to build all files generated by this target.
doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/build: doc/snippets/compile_TopicAliasing_block

.PHONY : doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/build

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/requires: doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/compile_TopicAliasing_block.cpp.o.requires

.PHONY : doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/requires

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/clean:
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_TopicAliasing_block.dir/cmake_clean.cmake
.PHONY : doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/clean

doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/depend:
	cd /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/doc/snippets /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets /home/chris/Documents/Repos/bicgstab_l/eigen_test_folder/build/doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/snippets/CMakeFiles/compile_TopicAliasing_block.dir/depend

