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
include src/CMakeFiles/gxlgedlib.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/gxlgedlib.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/gxlgedlib.dir/flags.make

src/CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.o: src/CMakeFiles/gxlgedlib.dir/flags.make
src/CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.o: ../src/env/ged_env.gxl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lucas/Documents/stage/gedlib/compression/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.o"
	cd /home/lucas/Documents/stage/gedlib/compression/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.o -c /home/lucas/Documents/stage/gedlib/src/env/ged_env.gxl.cpp

src/CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.i"
	cd /home/lucas/Documents/stage/gedlib/compression/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lucas/Documents/stage/gedlib/src/env/ged_env.gxl.cpp > CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.i

src/CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.s"
	cd /home/lucas/Documents/stage/gedlib/compression/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lucas/Documents/stage/gedlib/src/env/ged_env.gxl.cpp -o CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.s

# Object files for target gxlgedlib
gxlgedlib_OBJECTS = \
"CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.o"

# External object files for target gxlgedlib
gxlgedlib_EXTERNAL_OBJECTS =

../lib/libgxlgedlib.so: src/CMakeFiles/gxlgedlib.dir/env/ged_env.gxl.cpp.o
../lib/libgxlgedlib.so: src/CMakeFiles/gxlgedlib.dir/build.make
../lib/libgxlgedlib.so: src/CMakeFiles/gxlgedlib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lucas/Documents/stage/gedlib/compression/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../lib/libgxlgedlib.so"
	cd /home/lucas/Documents/stage/gedlib/compression/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gxlgedlib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/gxlgedlib.dir/build: ../lib/libgxlgedlib.so

.PHONY : src/CMakeFiles/gxlgedlib.dir/build

src/CMakeFiles/gxlgedlib.dir/clean:
	cd /home/lucas/Documents/stage/gedlib/compression/src && $(CMAKE_COMMAND) -P CMakeFiles/gxlgedlib.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/gxlgedlib.dir/clean

src/CMakeFiles/gxlgedlib.dir/depend:
	cd /home/lucas/Documents/stage/gedlib/compression && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lucas/Documents/stage/gedlib /home/lucas/Documents/stage/gedlib/src /home/lucas/Documents/stage/gedlib/compression /home/lucas/Documents/stage/gedlib/compression/src /home/lucas/Documents/stage/gedlib/compression/src/CMakeFiles/gxlgedlib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/gxlgedlib.dir/depend
