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
CMAKE_BINARY_DIR = /home/lucas/Documents/stage/gedlib/build

# Utility rule file for vldbj2020.

# Include the progress variables for this target.
include tests/vldbj2020/CMakeFiles/vldbj2020.dir/progress.make

tests/vldbj2020/CMakeFiles/vldbj2020: ../tests/vldbj2020/bin/vldbj_train_ring
tests/vldbj2020/CMakeFiles/vldbj2020: ../tests/vldbj2020/bin/vldbj_train_walks
tests/vldbj2020/CMakeFiles/vldbj2020: ../tests/vldbj2020/bin/vldbj_train_subgraph
tests/vldbj2020/CMakeFiles/vldbj2020: ../tests/vldbj2020/bin/vldbj_train_ml
tests/vldbj2020/CMakeFiles/vldbj2020: ../tests/vldbj2020/bin/vldbj_test_lsape_based_methods
tests/vldbj2020/CMakeFiles/vldbj2020: ../tests/vldbj2020/bin/vldbj_test_best_methods


vldbj2020: tests/vldbj2020/CMakeFiles/vldbj2020
vldbj2020: tests/vldbj2020/CMakeFiles/vldbj2020.dir/build.make

.PHONY : vldbj2020

# Rule to build all files generated by this target.
tests/vldbj2020/CMakeFiles/vldbj2020.dir/build: vldbj2020

.PHONY : tests/vldbj2020/CMakeFiles/vldbj2020.dir/build

tests/vldbj2020/CMakeFiles/vldbj2020.dir/clean:
	cd /home/lucas/Documents/stage/gedlib/build/tests/vldbj2020 && $(CMAKE_COMMAND) -P CMakeFiles/vldbj2020.dir/cmake_clean.cmake
.PHONY : tests/vldbj2020/CMakeFiles/vldbj2020.dir/clean

tests/vldbj2020/CMakeFiles/vldbj2020.dir/depend:
	cd /home/lucas/Documents/stage/gedlib/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Documents/stage/gedlib /home/lucas/Documents/stage/gedlib/tests/vldbj2020 /home/lucas/Documents/stage/gedlib/build /home/lucas/Documents/stage/gedlib/build/tests/vldbj2020 /home/lucas/Documents/stage/gedlib/build/tests/vldbj2020/CMakeFiles/vldbj2020.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/vldbj2020/CMakeFiles/vldbj2020.dir/depend

