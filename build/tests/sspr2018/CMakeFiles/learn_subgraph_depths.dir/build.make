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

# Include any dependencies generated for this target.
include tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/depend.make

# Include the progress variables for this target.
include tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/progress.make

# Include the compile flags for this target's objects.
include tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/flags.make

tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.o: tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/flags.make
tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.o: ../tests/sspr2018/src/learn_subgraph_depths.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucas/Documents/stage/gedlib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.o"
	cd /home/lucas/Documents/stage/gedlib/build/tests/sspr2018 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.o -c /home/lucas/Documents/stage/gedlib/tests/sspr2018/src/learn_subgraph_depths.cpp

tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.i"
	cd /home/lucas/Documents/stage/gedlib/build/tests/sspr2018 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lucas/Documents/stage/gedlib/tests/sspr2018/src/learn_subgraph_depths.cpp > CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.i

tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.s"
	cd /home/lucas/Documents/stage/gedlib/build/tests/sspr2018 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lucas/Documents/stage/gedlib/tests/sspr2018/src/learn_subgraph_depths.cpp -o CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.s

# Object files for target learn_subgraph_depths
learn_subgraph_depths_OBJECTS = \
"CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.o"

# External object files for target learn_subgraph_depths
learn_subgraph_depths_EXTERNAL_OBJECTS =

../tests/sspr2018/bin/learn_subgraph_depths: tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/src/learn_subgraph_depths.cpp.o
../tests/sspr2018/bin/learn_subgraph_depths: tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/build.make
../tests/sspr2018/bin/learn_subgraph_depths: ../lib/libgxlgedlib.so
../tests/sspr2018/bin/learn_subgraph_depths: tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucas/Documents/stage/gedlib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../tests/sspr2018/bin/learn_subgraph_depths"
	cd /home/lucas/Documents/stage/gedlib/build/tests/sspr2018 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/learn_subgraph_depths.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/build: ../tests/sspr2018/bin/learn_subgraph_depths

.PHONY : tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/build

tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/clean:
	cd /home/lucas/Documents/stage/gedlib/build/tests/sspr2018 && $(CMAKE_COMMAND) -P CMakeFiles/learn_subgraph_depths.dir/cmake_clean.cmake
.PHONY : tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/clean

tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/depend:
	cd /home/lucas/Documents/stage/gedlib/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Documents/stage/gedlib /home/lucas/Documents/stage/gedlib/tests/sspr2018 /home/lucas/Documents/stage/gedlib/build /home/lucas/Documents/stage/gedlib/build/tests/sspr2018 /home/lucas/Documents/stage/gedlib/build/tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/sspr2018/CMakeFiles/learn_subgraph_depths.dir/depend

