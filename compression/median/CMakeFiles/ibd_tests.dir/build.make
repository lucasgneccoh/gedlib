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
include median/CMakeFiles/ibd_tests.dir/depend.make

# Include the progress variables for this target.
include median/CMakeFiles/ibd_tests.dir/progress.make

# Include the compile flags for this target's objects.
include median/CMakeFiles/ibd_tests.dir/flags.make

median/CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.o: median/CMakeFiles/ibd_tests.dir/flags.make
median/CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.o: ../median/tests/ibd_tests.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucas/Documents/stage/gedlib/compression/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object median/CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.o"
	cd /home/lucas/Documents/stage/gedlib/compression/median && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.o -c /home/lucas/Documents/stage/gedlib/median/tests/ibd_tests.cpp

median/CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.i"
	cd /home/lucas/Documents/stage/gedlib/compression/median && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lucas/Documents/stage/gedlib/median/tests/ibd_tests.cpp > CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.i

median/CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.s"
	cd /home/lucas/Documents/stage/gedlib/compression/median && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lucas/Documents/stage/gedlib/median/tests/ibd_tests.cpp -o CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.s

# Object files for target ibd_tests
ibd_tests_OBJECTS = \
"CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.o"

# External object files for target ibd_tests
ibd_tests_EXTERNAL_OBJECTS =

../median/bin/ibd_tests: median/CMakeFiles/ibd_tests.dir/tests/ibd_tests.cpp.o
../median/bin/ibd_tests: median/CMakeFiles/ibd_tests.dir/build.make
../median/bin/ibd_tests: ../lib/libgxlgedlib.so
../median/bin/ibd_tests: median/CMakeFiles/ibd_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucas/Documents/stage/gedlib/compression/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../median/bin/ibd_tests"
	cd /home/lucas/Documents/stage/gedlib/compression/median && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ibd_tests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
median/CMakeFiles/ibd_tests.dir/build: ../median/bin/ibd_tests

.PHONY : median/CMakeFiles/ibd_tests.dir/build

median/CMakeFiles/ibd_tests.dir/clean:
	cd /home/lucas/Documents/stage/gedlib/compression/median && $(CMAKE_COMMAND) -P CMakeFiles/ibd_tests.dir/cmake_clean.cmake
.PHONY : median/CMakeFiles/ibd_tests.dir/clean

median/CMakeFiles/ibd_tests.dir/depend:
	cd /home/lucas/Documents/stage/gedlib/compression && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Documents/stage/gedlib /home/lucas/Documents/stage/gedlib/median /home/lucas/Documents/stage/gedlib/compression /home/lucas/Documents/stage/gedlib/compression/median /home/lucas/Documents/stage/gedlib/compression/median/CMakeFiles/ibd_tests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : median/CMakeFiles/ibd_tests.dir/depend
