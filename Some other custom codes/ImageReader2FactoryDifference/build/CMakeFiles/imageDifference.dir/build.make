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
CMAKE_SOURCE_DIR = /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/build

# Include any dependencies generated for this target.
include CMakeFiles/imageDifference.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/imageDifference.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/imageDifference.dir/flags.make

CMakeFiles/imageDifference.dir/imageDifference.cxx.o: CMakeFiles/imageDifference.dir/flags.make
CMakeFiles/imageDifference.dir/imageDifference.cxx.o: ../imageDifference.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/imageDifference.dir/imageDifference.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/imageDifference.dir/imageDifference.cxx.o -c /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/imageDifference.cxx

CMakeFiles/imageDifference.dir/imageDifference.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imageDifference.dir/imageDifference.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/imageDifference.cxx > CMakeFiles/imageDifference.dir/imageDifference.cxx.i

CMakeFiles/imageDifference.dir/imageDifference.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imageDifference.dir/imageDifference.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/imageDifference.cxx -o CMakeFiles/imageDifference.dir/imageDifference.cxx.s

CMakeFiles/imageDifference.dir/imageDifference.cxx.o.requires:
.PHONY : CMakeFiles/imageDifference.dir/imageDifference.cxx.o.requires

CMakeFiles/imageDifference.dir/imageDifference.cxx.o.provides: CMakeFiles/imageDifference.dir/imageDifference.cxx.o.requires
	$(MAKE) -f CMakeFiles/imageDifference.dir/build.make CMakeFiles/imageDifference.dir/imageDifference.cxx.o.provides.build
.PHONY : CMakeFiles/imageDifference.dir/imageDifference.cxx.o.provides

CMakeFiles/imageDifference.dir/imageDifference.cxx.o.provides.build: CMakeFiles/imageDifference.dir/imageDifference.cxx.o

# Object files for target imageDifference
imageDifference_OBJECTS = \
"CMakeFiles/imageDifference.dir/imageDifference.cxx.o"

# External object files for target imageDifference
imageDifference_EXTERNAL_OBJECTS =

imageDifference.app/Contents/MacOS/imageDifference: CMakeFiles/imageDifference.dir/imageDifference.cxx.o
imageDifference.app/Contents/MacOS/imageDifference: CMakeFiles/imageDifference.dir/build.make
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkDomainsChemistry-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersFlowPaths-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersGeneric-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersHyperTree-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersParallelImaging-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersProgrammable-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersSelection-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersSMP-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersTexture-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersVerdict-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkverdict-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkGeovisCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkproj4-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingMath-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingMorphological-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingStatistics-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingStencil-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkInteractionImage-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOAMR-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOEnSight-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOExodus-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOExport-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingGL2PS-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingContextOpenGL-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOImport-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOInfovis-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtklibxml2-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOLSDyna-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOMINC-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOMovie-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkoggtheora-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOParallel-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOParallelXML-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOPLY-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOVideo-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingFreeTypeOpenGL-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingImage-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingLIC-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingLOD-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingVolumeOpenGL-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkTestingGenericBridge-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkTestingIOSQL-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkTestingRendering-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkViewsContext2D-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkViewsInfovis-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersAMR-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkgl2ps-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkexoIIc-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersParallel-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIONetCDF-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkNetCDF_cxx-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkNetCDF-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkhdf5_hl-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkhdf5-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkParallelCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOXML-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOGeometry-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkjsoncpp-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOXMLParser-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkexpat-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOLegacy-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingOpenGL-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOSQL-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtksqlite-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkChartsCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingContext2D-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersImaging-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkInfovisLayout-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkInfovisCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkViewsCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkInteractionWidgets-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersHybrid-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingGeneral-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingSources-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersModeling-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingHybrid-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOImage-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkDICOMParser-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkIOCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkmetaio-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkpng-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtktiff-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkjpeg-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkInteractionStyle-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingAnnotation-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingColor-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingVolume-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingLabel-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingFreeType-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkRenderingCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonColor-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersExtraction-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersStatistics-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkalglib-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingFourier-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkImagingCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersGeometry-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersSources-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersGeneral-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkFiltersCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonExecutionModel-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonComputationalGeometry-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonDataModel-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonMisc-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonTransforms-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonMath-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonSystem-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkCommonCore-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtksys-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkftgl-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkfreetype-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: /usr/local/lib/libvtkzlib-6.2.1.dylib
imageDifference.app/Contents/MacOS/imageDifference: CMakeFiles/imageDifference.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable imageDifference.app/Contents/MacOS/imageDifference"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/imageDifference.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/imageDifference.dir/build: imageDifference.app/Contents/MacOS/imageDifference
.PHONY : CMakeFiles/imageDifference.dir/build

CMakeFiles/imageDifference.dir/requires: CMakeFiles/imageDifference.dir/imageDifference.cxx.o.requires
.PHONY : CMakeFiles/imageDifference.dir/requires

CMakeFiles/imageDifference.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/imageDifference.dir/cmake_clean.cmake
.PHONY : CMakeFiles/imageDifference.dir/clean

CMakeFiles/imageDifference.dir/depend:
	cd /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/build /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/build /Users/ashraymalhotra/Downloads/ImageReader2FactoryDifference/build/CMakeFiles/imageDifference.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/imageDifference.dir/depend

