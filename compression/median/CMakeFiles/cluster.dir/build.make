# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/lucas/Documents/stage/gedlib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lucas/Documents/stage/gedlib/compression

# Utility rule file for cluster.

# Include the progress variables for this target.
include median/CMakeFiles/cluster.dir/progress.make

median/CMakeFiles/cluster: ../median/bin/cluster_letter
median/CMakeFiles/cluster: ../median/bin/clustering_tests
median/CMakeFiles/cluster: ../median/bin/classification_tests


cluster: median/CMakeFiles/cluster
cluster: median/CMakeFiles/cluster.dir/build.make

.PHONY : cluster

# Rule to build all files generated by this target.
median/CMakeFiles/cluster.dir/build: cluster

.PHONY : median/CMakeFiles/cluster.dir/build

median/CMakeFiles/cluster.dir/clean:
	cd /home/lucas/Documents/stage/gedlib/compression/median && $(CMAKE_COMMAND) -P CMakeFiles/cluster.dir/cmake_clean.cmake
.PHONY : median/CMakeFiles/cluster.dir/clean

median/CMakeFiles/cluster.dir/depend:
	cd /home/lucas/Documents/stage/gedlib/compression && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Documents/stage/gedlib /home/lucas/Documents/stage/gedlib/median /home/lucas/Documents/stage/gedlib/compression /home/lucas/Documents/stage/gedlib/compression/median /home/lucas/Documents/stage/gedlib/compression/median/CMakeFiles/cluster.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : median/CMakeFiles/cluster.dir/depend

