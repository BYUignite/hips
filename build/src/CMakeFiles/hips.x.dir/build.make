# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.23.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.23.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/me/codes/reactHips/hips

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/me/codes/reactHips/hips/build

# Include any dependencies generated for this target.
include src/CMakeFiles/hips.x.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/hips.x.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/hips.x.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/hips.x.dir/flags.make

src/CMakeFiles/hips.x.dir/driver.cc.o: src/CMakeFiles/hips.x.dir/flags.make
src/CMakeFiles/hips.x.dir/driver.cc.o: ../src/driver.cc
src/CMakeFiles/hips.x.dir/driver.cc.o: src/CMakeFiles/hips.x.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/me/codes/reactHips/hips/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/hips.x.dir/driver.cc.o"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/hips.x.dir/driver.cc.o -MF CMakeFiles/hips.x.dir/driver.cc.o.d -o CMakeFiles/hips.x.dir/driver.cc.o -c /Users/me/codes/reactHips/hips/src/driver.cc

src/CMakeFiles/hips.x.dir/driver.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hips.x.dir/driver.cc.i"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/me/codes/reactHips/hips/src/driver.cc > CMakeFiles/hips.x.dir/driver.cc.i

src/CMakeFiles/hips.x.dir/driver.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hips.x.dir/driver.cc.s"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/me/codes/reactHips/hips/src/driver.cc -o CMakeFiles/hips.x.dir/driver.cc.s

src/CMakeFiles/hips.x.dir/cvodeDriver.cc.o: src/CMakeFiles/hips.x.dir/flags.make
src/CMakeFiles/hips.x.dir/cvodeDriver.cc.o: ../src/cvodeDriver.cc
src/CMakeFiles/hips.x.dir/cvodeDriver.cc.o: src/CMakeFiles/hips.x.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/me/codes/reactHips/hips/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/hips.x.dir/cvodeDriver.cc.o"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/hips.x.dir/cvodeDriver.cc.o -MF CMakeFiles/hips.x.dir/cvodeDriver.cc.o.d -o CMakeFiles/hips.x.dir/cvodeDriver.cc.o -c /Users/me/codes/reactHips/hips/src/cvodeDriver.cc

src/CMakeFiles/hips.x.dir/cvodeDriver.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hips.x.dir/cvodeDriver.cc.i"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/me/codes/reactHips/hips/src/cvodeDriver.cc > CMakeFiles/hips.x.dir/cvodeDriver.cc.i

src/CMakeFiles/hips.x.dir/cvodeDriver.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hips.x.dir/cvodeDriver.cc.s"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/me/codes/reactHips/hips/src/cvodeDriver.cc -o CMakeFiles/hips.x.dir/cvodeDriver.cc.s

src/CMakeFiles/hips.x.dir/hips.cc.o: src/CMakeFiles/hips.x.dir/flags.make
src/CMakeFiles/hips.x.dir/hips.cc.o: ../src/hips.cc
src/CMakeFiles/hips.x.dir/hips.cc.o: src/CMakeFiles/hips.x.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/me/codes/reactHips/hips/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/hips.x.dir/hips.cc.o"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/hips.x.dir/hips.cc.o -MF CMakeFiles/hips.x.dir/hips.cc.o.d -o CMakeFiles/hips.x.dir/hips.cc.o -c /Users/me/codes/reactHips/hips/src/hips.cc

src/CMakeFiles/hips.x.dir/hips.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hips.x.dir/hips.cc.i"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/me/codes/reactHips/hips/src/hips.cc > CMakeFiles/hips.x.dir/hips.cc.i

src/CMakeFiles/hips.x.dir/hips.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hips.x.dir/hips.cc.s"
	cd /Users/me/codes/reactHips/hips/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/me/codes/reactHips/hips/src/hips.cc -o CMakeFiles/hips.x.dir/hips.cc.s

# Object files for target hips.x
hips_x_OBJECTS = \
"CMakeFiles/hips.x.dir/driver.cc.o" \
"CMakeFiles/hips.x.dir/cvodeDriver.cc.o" \
"CMakeFiles/hips.x.dir/hips.cc.o"

# External object files for target hips.x
hips_x_EXTERNAL_OBJECTS =

src/hips.x: src/CMakeFiles/hips.x.dir/driver.cc.o
src/hips.x: src/CMakeFiles/hips.x.dir/cvodeDriver.cc.o
src/hips.x: src/CMakeFiles/hips.x.dir/hips.cc.o
src/hips.x: src/CMakeFiles/hips.x.dir/build.make
src/hips.x: /usr/local/lib/libcantera.a
src/hips.x: /usr/local/lib/libsundials_cvodes.dylib
src/hips.x: /usr/local/opt/python@3.10/Frameworks/Python.framework/Versions/3.10/lib/libpython3.10.dylib
src/hips.x: src/CMakeFiles/hips.x.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/me/codes/reactHips/hips/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable hips.x"
	cd /Users/me/codes/reactHips/hips/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hips.x.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/hips.x.dir/build: src/hips.x
.PHONY : src/CMakeFiles/hips.x.dir/build

src/CMakeFiles/hips.x.dir/clean:
	cd /Users/me/codes/reactHips/hips/build/src && $(CMAKE_COMMAND) -P CMakeFiles/hips.x.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/hips.x.dir/clean

src/CMakeFiles/hips.x.dir/depend:
	cd /Users/me/codes/reactHips/hips/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/me/codes/reactHips/hips /Users/me/codes/reactHips/hips/src /Users/me/codes/reactHips/hips/build /Users/me/codes/reactHips/hips/build/src /Users/me/codes/reactHips/hips/build/src/CMakeFiles/hips.x.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/hips.x.dir/depend

