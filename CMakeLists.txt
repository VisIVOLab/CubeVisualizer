#cmake -G Xcode  -DCFITSIO_INC_PATH=/opt/cfitsio-3.49/include -DCFITSIO_LIB_PATH=/opt/cfitsio-3.49/lib -DBOOST_INC_PATH=/opt/boost_1_75_0 .. 


cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(FitsCubeDemo LANGUAGES CXX)

# VTK 9.0.3
find_package(VTK COMPONENTS 
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonTransforms
  vtkFiltersCore
  vtkFiltersModeling
  vtkFiltersSources
  vtkIOImage
  vtkImagingHybrid
  vtkInteractionStyle
  vtkInteractionWidgets
  vtkRenderingAnnotation
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkRenderingVolume
  vtkRenderingVolumeOpenGL2
  QUIET
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Skipping FitsCubeDemo: ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()

# CFITSIO 3.49
find_library (CFITSIO_LIB NAMES cfitsio PATHS ${CFITSIO_LIB_PATH})
include_directories(${CFITSIO_INC_PATH})

# BOOST 1.75.0
include_directories(${BOOST_INC_PATH})

# OpenMP
find_library (OPENMP_LIB NAMES omp PATHS ${OPENMP_LIB_PATH})
include_directories(${OPENMP_INC_PATH})
set(CMAKE_CXX_FLAGS "-Xpreprocessor -fopenmp")

message (STATUS "VTK_VERSION: ${VTK_VERSION}")
message (STATUS "VTK_LIBRARIES: ${VTK_LIBRARIES}")
message (STATUS "CFITSIO_LIB: ${CFITSIO_LIB}")
message (STATUS "CFITSIO_INC_PATH: ${CFITSIO_INC_PATH}")
message (STATUS "BOOST_INC_PATH: ${BOOST_INC_PATH}")
message (STATUS "OpenMP_LIB: ${OPENMP_LIB}")
message (STATUS "OpenMP_INC_PATH: ${OPENMP_INC_PATH}")

# include all components
add_executable(FitsCubeDemo MACOSX_BUNDLE FitsCubeDemo.cxx vtkfitsreader.cpp vtkfitsreader.h)
target_link_libraries(FitsCubeDemo PRIVATE ${CFITSIO_LIB} ${VTK_LIBRARIES})
target_link_libraries(FitsCubeDemo PRIVATE ${OPENMP_LIB})

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS FitsCubeDemo
  MODULES ${VTK_LIBRARIES}
)
