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
include tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/depend.make

# Include the progress variables for this target.
include tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/progress.make

# Include the compile flags for this target's objects.
include tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/flags.make

tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.o: tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/flags.make
tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.o: ../tests/vldbj2020/src/train_ring.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucas/Documents/stage/gedlib/compression/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.o"
	cd /home/lucas/Documents/stage/gedlib/compression/tests/vldbj2020 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.o -c /home/lucas/Documents/stage/gedlib/tests/vldbj2020/src/train_ring.cpp

tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.i"
	cd /home/lucas/Documents/stage/gedlib/compression/tests/vldbj2020 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lucas/Documents/stage/gedlib/tests/vldbj2020/src/train_ring.cpp > CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.i

tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.s"
	cd /home/lucas/Documents/stage/gedlib/compression/tests/vldbj2020 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lucas/Documents/stage/gedlib/tests/vldbj2020/src/train_ring.cpp -o CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.s

# Object files for target vldbj_train_ring
vldbj_train_ring_OBJECTS = \
"CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.o"

# External object files for target vldbj_train_ring
vldbj_train_ring_EXTERNAL_OBJECTS =

../tests/vldbj2020/bin/vldbj_train_ring: tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/src/train_ring.cpp.o
../tests/vldbj2020/bin/vldbj_train_ring: tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/build.make
../tests/vldbj2020/bin/vldbj_train_ring: ../lib/libgxlgedlib.so
../tests/vldbj2020/bin/vldbj_train_ring: tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucas/Documents/stage/gedlib/compression/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../tests/vldbj2020/bin/vldbj_train_ring"
	cd /home/lucas/Documents/stage/gedlib/compression/tests/vldbj2020 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vldbj_train_ring.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/build: ../tests/vldbj2020/bin/vldbj_train_ring

.PHONY : tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/build

tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/clean:
	cd /home/lucas/Documents/stage/gedlib/compression/tests/vldbj2020 && $(CMAKE_COMMAND) -P CMakeFiles/vldbj_train_ring.dir/cmake_clean.cmake
.PHONY : tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/clean

tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/depend:
	cd /home/lucas/Documents/stage/gedlib/compression && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Documents/stage/gedlib /home/lucas/Documents/stage/gedlib/tests/vldbj2020 /home/lucas/Documents/stage/gedlib/compression /home/lucas/Documents/stage/gedlib/compression/tests/vldbj2020 /home/lucas/Documents/stage/gedlib/compression/tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/vldbj2020/CMakeFiles/vldbj_train_ring.dir/depend

