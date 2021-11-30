#cmake -G Xcode  -DCFITSIO_INC_PATH=/opt/cfitsio-3.49/include -DCFITSIO_LIB_PATH=/opt/cfitsio-3.49/lib -DBOOST_INC_PATH=/opt/boost_1_75_0 .. 


cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(FitsCubeDemo)

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

find_library (CFITSIO_LIB NAMES cfitsio PATHS ${CFITSIO_LIB_PATH})
include_directories(${CFITSIO_INC_PATH})
include_directories(${BOOST_INC_PATH})



if (NOT VTK_FOUND)
  message("Skipping FitsCubeDemo: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
message (STATUS "CFITSIO_LIB: ${CFITSIO_LIB}")
message (STATUS "CFITSIO_INC_PATH: ${CFITSIO_INC_PATH}")
message (STATUS "BOOST_INC_PATH: ${BOOST_INC_PATH}")
message (STATUS "VTK_LIBRARIES: ${VTK_LIBRARIES}")

# include all components
add_executable(FitsCubeDemo MACOSX_BUNDLE FitsCubeDemo.cxx vtkfitsreader.cpp)
target_link_libraries(FitsCubeDemo PRIVATE ${CFITSIO_LIB} ${VTK_LIBRARIES})
# vtk_module_autoinit is needed
vtk_module_autoinit(
	TARGETS FitsCubeDemo
    MODULES ${VTK_LIBRARIES}
)
