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
CMAKE_SOURCE_DIR = "/home/andruon/Desktop/Elem finis/Integrate"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/andruon/Desktop/Elem finis/Integrate/build"

# Include any dependencies generated for this target.
include deps/bov/CMakeFiles/koch_snowflake.dir/depend.make

# Include the progress variables for this target.
include deps/bov/CMakeFiles/koch_snowflake.dir/progress.make

# Include the compile flags for this target's objects.
include deps/bov/CMakeFiles/koch_snowflake.dir/flags.make

deps/bov/CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.o: deps/bov/CMakeFiles/koch_snowflake.dir/flags.make
deps/bov/CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.o: ../deps/bov/examples/koch_snowflake.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/andruon/Desktop/Elem finis/Integrate/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object deps/bov/CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.o"
	cd "/home/andruon/Desktop/Elem finis/Integrate/build/deps/bov" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.o   -c "/home/andruon/Desktop/Elem finis/Integrate/deps/bov/examples/koch_snowflake.c"

deps/bov/CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.i"
	cd "/home/andruon/Desktop/Elem finis/Integrate/build/deps/bov" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/andruon/Desktop/Elem finis/Integrate/deps/bov/examples/koch_snowflake.c" > CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.i

deps/bov/CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.s"
	cd "/home/andruon/Desktop/Elem finis/Integrate/build/deps/bov" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/andruon/Desktop/Elem finis/Integrate/deps/bov/examples/koch_snowflake.c" -o CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.s

# Object files for target koch_snowflake
koch_snowflake_OBJECTS = \
"CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.o"

# External object files for target koch_snowflake
koch_snowflake_EXTERNAL_OBJECTS =

deps/bov/examples/koch_snowflake: deps/bov/CMakeFiles/koch_snowflake.dir/examples/koch_snowflake.c.o
deps/bov/examples/koch_snowflake: deps/bov/CMakeFiles/koch_snowflake.dir/build.make
deps/bov/examples/koch_snowflake: deps/bov/lib/libbov.a
deps/bov/examples/koch_snowflake: deps/bov/deps/glad/libglad.a
deps/bov/examples/koch_snowflake: deps/bov/deps/glfw/src/libglfw3.a
deps/bov/examples/koch_snowflake: /usr/lib/x86_64-linux-gnu/librt.so
deps/bov/examples/koch_snowflake: /usr/lib/x86_64-linux-gnu/libm.so
deps/bov/examples/koch_snowflake: /usr/lib/x86_64-linux-gnu/libX11.so
deps/bov/examples/koch_snowflake: deps/bov/CMakeFiles/koch_snowflake.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/andruon/Desktop/Elem finis/Integrate/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable examples/koch_snowflake"
	cd "/home/andruon/Desktop/Elem finis/Integrate/build/deps/bov" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/koch_snowflake.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
deps/bov/CMakeFiles/koch_snowflake.dir/build: deps/bov/examples/koch_snowflake

.PHONY : deps/bov/CMakeFiles/koch_snowflake.dir/build

deps/bov/CMakeFiles/koch_snowflake.dir/clean:
	cd "/home/andruon/Desktop/Elem finis/Integrate/build/deps/bov" && $(CMAKE_COMMAND) -P CMakeFiles/koch_snowflake.dir/cmake_clean.cmake
.PHONY : deps/bov/CMakeFiles/koch_snowflake.dir/clean

deps/bov/CMakeFiles/koch_snowflake.dir/depend:
	cd "/home/andruon/Desktop/Elem finis/Integrate/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/andruon/Desktop/Elem finis/Integrate" "/home/andruon/Desktop/Elem finis/Integrate/deps/bov" "/home/andruon/Desktop/Elem finis/Integrate/build" "/home/andruon/Desktop/Elem finis/Integrate/build/deps/bov" "/home/andruon/Desktop/Elem finis/Integrate/build/deps/bov/CMakeFiles/koch_snowflake.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : deps/bov/CMakeFiles/koch_snowflake.dir/depend

