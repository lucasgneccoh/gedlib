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

# Include any dependencies generated for this target.
include tests/tkde2019/CMakeFiles/train_subgraph.dir/depend.make

# Include the progress variables for this target.
include tests/tkde2019/CMakeFiles/train_subgraph.dir/progress.make

# Include the compile flags for this target's objects.
include tests/tkde2019/CMakeFiles/train_subgraph.dir/flags.make

tests/tkde2019/CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.o: tests/tkde2019/CMakeFiles/train_subgraph.dir/flags.make
tests/tkde2019/CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.o: ../tests/tkde2019/src/train_subgraph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucas/Documents/stage/gedlib/compression/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/tkde2019/CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.o"
	cd /home/lucas/Documents/stage/gedlib/compression/tests/tkde2019 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.o -c /home/lucas/Documents/stage/gedlib/tests/tkde2019/src/train_subgraph.cpp

tests/tkde2019/CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.i"
	cd /home/lucas/Documents/stage/gedlib/compression/tests/tkde2019 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lucas/Documents/stage/gedlib/tests/tkde2019/src/train_subgraph.cpp > CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.i

tests/tkde2019/CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.s"
	cd /home/lucas/Documents/stage/gedlib/compression/tests/tkde2019 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lucas/Documents/stage/gedlib/tests/tkde2019/src/train_subgraph.cpp -o CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.s

# Object files for target train_subgraph
train_subgraph_OBJECTS = \
"CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.o"

# External object files for target train_subgraph
train_subgraph_EXTERNAL_OBJECTS =

../tests/tkde2019/bin/train_subgraph: tests/tkde2019/CMakeFiles/train_subgraph.dir/src/train_subgraph.cpp.o
../tests/tkde2019/bin/train_subgraph: tests/tkde2019/CMakeFiles/train_subgraph.dir/build.make
../tests/tkde2019/bin/train_subgraph: ../lib/libgxlgedlib.so
../tests/tkde2019/bin/train_subgraph: tests/tkde2019/CMakeFiles/train_subgraph.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucas/Documents/stage/gedlib/compression/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../tests/tkde2019/bin/train_subgraph"
	cd /home/lucas/Documents/stage/gedlib/compression/tests/tkde2019 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/train_subgraph.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/tkde2019/CMakeFiles/train_subgraph.dir/build: ../tests/tkde2019/bin/train_subgraph

.PHONY : tests/tkde2019/CMakeFiles/train_subgraph.dir/build

tests/tkde2019/CMakeFiles/train_subgraph.dir/clean:
	cd /home/lucas/Documents/stage/gedlib/compression/tests/tkde2019 && $(CMAKE_COMMAND) -P CMakeFiles/train_subgraph.dir/cmake_clean.cmake
.PHONY : tests/tkde2019/CMakeFiles/train_subgraph.dir/clean

tests/tkde2019/CMakeFiles/train_subgraph.dir/depend:
	cd /home/lucas/Documents/stage/gedlib/compression && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Documents/stage/gedlib /home/lucas/Documents/stage/gedlib/tests/tkde2019 /home/lucas/Documents/stage/gedlib/compression /home/lucas/Documents/stage/gedlib/compression/tests/tkde2019 /home/lucas/Documents/stage/gedlib/compression/tests/tkde2019/CMakeFiles/train_subgraph.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/tkde2019/CMakeFiles/train_subgraph.dir/depend

