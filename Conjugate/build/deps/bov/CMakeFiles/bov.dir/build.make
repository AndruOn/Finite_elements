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
CMAKE_SOURCE_DIR = "/home/andruon/Desktop/Elem finis/Conjugate"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/andruon/Desktop/Elem finis/Conjugate/build"

# Include any dependencies generated for this target.
include deps/bov/CMakeFiles/bov.dir/depend.make

# Include the progress variables for this target.
include deps/bov/CMakeFiles/bov.dir/progress.make

# Include the compile flags for this target's objects.
include deps/bov/CMakeFiles/bov.dir/flags.make

deps/bov/CMakeFiles/bov.dir/src/BOV.c.o: deps/bov/CMakeFiles/bov.dir/flags.make
deps/bov/CMakeFiles/bov.dir/src/BOV.c.o: ../deps/bov/src/BOV.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/andruon/Desktop/Elem finis/Conjugate/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object deps/bov/CMakeFiles/bov.dir/src/BOV.c.o"
	cd "/home/andruon/Desktop/Elem finis/Conjugate/build/deps/bov" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/bov.dir/src/BOV.c.o   -c "/home/andruon/Desktop/Elem finis/Conjugate/deps/bov/src/BOV.c"

deps/bov/CMakeFiles/bov.dir/src/BOV.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/bov.dir/src/BOV.c.i"
	cd "/home/andruon/Desktop/Elem finis/Conjugate/build/deps/bov" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/andruon/Desktop/Elem finis/Conjugate/deps/bov/src/BOV.c" > CMakeFiles/bov.dir/src/BOV.c.i

deps/bov/CMakeFiles/bov.dir/src/BOV.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/bov.dir/src/BOV.c.s"
	cd "/home/andruon/Desktop/Elem finis/Conjugate/build/deps/bov" && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/andruon/Desktop/Elem finis/Conjugate/deps/bov/src/BOV.c" -o CMakeFiles/bov.dir/src/BOV.c.s

# Object files for target bov
bov_OBJECTS = \
"CMakeFiles/bov.dir/src/BOV.c.o"

# External object files for target bov
bov_EXTERNAL_OBJECTS =

deps/bov/lib/libbov.a: deps/bov/CMakeFiles/bov.dir/src/BOV.c.o
deps/bov/lib/libbov.a: deps/bov/CMakeFiles/bov.dir/build.make
deps/bov/lib/libbov.a: deps/bov/CMakeFiles/bov.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/andruon/Desktop/Elem finis/Conjugate/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library lib/libbov.a"
	cd "/home/andruon/Desktop/Elem finis/Conjugate/build/deps/bov" && $(CMAKE_COMMAND) -P CMakeFiles/bov.dir/cmake_clean_target.cmake
	cd "/home/andruon/Desktop/Elem finis/Conjugate/build/deps/bov" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bov.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
deps/bov/CMakeFiles/bov.dir/build: deps/bov/lib/libbov.a

.PHONY : deps/bov/CMakeFiles/bov.dir/build

deps/bov/CMakeFiles/bov.dir/clean:
	cd "/home/andruon/Desktop/Elem finis/Conjugate/build/deps/bov" && $(CMAKE_COMMAND) -P CMakeFiles/bov.dir/cmake_clean.cmake
.PHONY : deps/bov/CMakeFiles/bov.dir/clean

deps/bov/CMakeFiles/bov.dir/depend:
	cd "/home/andruon/Desktop/Elem finis/Conjugate/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/andruon/Desktop/Elem finis/Conjugate" "/home/andruon/Desktop/Elem finis/Conjugate/deps/bov" "/home/andruon/Desktop/Elem finis/Conjugate/build" "/home/andruon/Desktop/Elem finis/Conjugate/build/deps/bov" "/home/andruon/Desktop/Elem finis/Conjugate/build/deps/bov/CMakeFiles/bov.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : deps/bov/CMakeFiles/bov.dir/depend

