# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.1

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.1.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.1.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ashraymalhotra/Downloads/ImageReader2Factory2images

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/build

# Include any dependencies generated for this target.
include CMakeFiles/2ImageReader.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/2ImageReader.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/2ImageReader.dir/flags.make

CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o: CMakeFiles/2ImageReader.dir/flags.make
CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o: ../2ImageReader.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o -c /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/2ImageReader.cxx

CMakeFiles/2ImageReader.dir/2ImageReader.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2ImageReader.dir/2ImageReader.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/2ImageReader.cxx > CMakeFiles/2ImageReader.dir/2ImageReader.cxx.i

CMakeFiles/2ImageReader.dir/2ImageReader.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2ImageReader.dir/2ImageReader.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/2ImageReader.cxx -o CMakeFiles/2ImageReader.dir/2ImageReader.cxx.s

CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o.requires:
.PHONY : CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o.requires

CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o.provides: CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o.requires
	$(MAKE) -f CMakeFiles/2ImageReader.dir/build.make CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o.provides.build
.PHONY : CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o.provides

CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o.provides.build: CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o

# Object files for target 2ImageReader
2ImageReader_OBJECTS = \
"CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o"

# External object files for target 2ImageReader
2ImageReader_EXTERNAL_OBJECTS =

2ImageReader.app/Contents/MacOS/2ImageReader: CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o
2ImageReader.app/Contents/MacOS/2ImageReader: CMakeFiles/2ImageReader.dir/build.make
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkDomainsChemistry-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersFlowPaths-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersGeneric-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersHyperTree-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersParallelImaging-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersProgrammable-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersSelection-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersSMP-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersTexture-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersVerdict-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkverdict-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkGeovisCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkproj4-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingMath-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingMorphological-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingStatistics-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingStencil-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkInteractionImage-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOAMR-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOEnSight-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOExodus-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOExport-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingGL2PS-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingContextOpenGL-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOImport-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOInfovis-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtklibxml2-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOLSDyna-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOMINC-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOMovie-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkoggtheora-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOParallel-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOParallelXML-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOPLY-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOVideo-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingFreeTypeOpenGL-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingImage-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingLIC-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingLOD-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingVolumeOpenGL-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkTestingGenericBridge-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkTestingIOSQL-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkTestingRendering-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkViewsContext2D-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkViewsInfovis-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersAMR-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkgl2ps-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkexoIIc-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersParallel-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIONetCDF-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkNetCDF_cxx-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkNetCDF-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkhdf5_hl-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkhdf5-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkParallelCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOXML-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOGeometry-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkjsoncpp-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOXMLParser-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkexpat-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOLegacy-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingOpenGL-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOSQL-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtksqlite-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkChartsCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingContext2D-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersImaging-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkInfovisLayout-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkInfovisCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkViewsCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkInteractionWidgets-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersHybrid-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingGeneral-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingSources-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersModeling-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingHybrid-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOImage-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkDICOMParser-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkIOCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkmetaio-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkpng-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtktiff-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkjpeg-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkInteractionStyle-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingAnnotation-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingColor-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingVolume-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingLabel-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingFreeType-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkRenderingCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonColor-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersExtraction-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersStatistics-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkalglib-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingFourier-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkImagingCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersGeometry-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersSources-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersGeneral-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkFiltersCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonExecutionModel-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonComputationalGeometry-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonDataModel-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonMisc-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonTransforms-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonMath-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonSystem-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkCommonCore-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtksys-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkftgl-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkfreetype-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: /usr/local/lib/libvtkzlib-6.2.1.dylib
2ImageReader.app/Contents/MacOS/2ImageReader: CMakeFiles/2ImageReader.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable 2ImageReader.app/Contents/MacOS/2ImageReader"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/2ImageReader.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/2ImageReader.dir/build: 2ImageReader.app/Contents/MacOS/2ImageReader
.PHONY : CMakeFiles/2ImageReader.dir/build

CMakeFiles/2ImageReader.dir/requires: CMakeFiles/2ImageReader.dir/2ImageReader.cxx.o.requires
.PHONY : CMakeFiles/2ImageReader.dir/requires

CMakeFiles/2ImageReader.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/2ImageReader.dir/cmake_clean.cmake
.PHONY : CMakeFiles/2ImageReader.dir/clean

CMakeFiles/2ImageReader.dir/depend:
	cd /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ashraymalhotra/Downloads/ImageReader2Factory2images /Users/ashraymalhotra/Downloads/ImageReader2Factory2images /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/build /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/build /Users/ashraymalhotra/Downloads/ImageReader2Factory2images/build/CMakeFiles/2ImageReader.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/2ImageReader.dir/depend

