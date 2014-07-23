cmake_minimum_required( VERSION 2.6 )

# create a basic project
project( pyPolyCSG )

# add the polyhcsg subdirectory which contains the polyhedron class
add_subdirectory(libpolyhcsg)

set( BOOLEAN_SOURCES
  source/pypolyhedron.cpp
)

set( BOOLEAN_HEADERS
  include/pypolyhedron.h
)

find_package( PythonLibs REQUIRED )
IF( PYTHONLIBS_FOUND )
  set( INCLUDE_DIRS ${INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} )
  set( LIBS ${LIBS} ${PYTHON_LIBRARIES} )
ENDIF( PYTHONLIBS_FOUND )

find_package( Boost COMPONENTS python thread REQUIRED )
IF( Boost_FOUND )
  SET( INCLUDE_DIRS ${INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} )
  SET( LIBS ${LIBS} ${Boost_LIBRARIES} )
ENDIF( Boost_FOUND )

IF( PYTHONLIBS_FOUND AND Boost_FOUND )
  set( BOOLEAN_SOURCES ${BOOLEAN_SOURCES} source/python_wrapper.cpp )
ENDIF( PYTHONLIBS_FOUND AND Boost_FOUND )

include_directories( include libpolyhcsg/include ${INCLUDE_DIRS} )
add_library( pyPolyCSG SHARED ${BOOLEAN_SOURCES} ${BOOLEAN_HEADERS} )
target_link_libraries( pyPolyCSG polyhcsg )
IF( APPLE )
  # OS-X specific library naming
  #set_target_properties( pyPolyCSG PROPERTIES PREFIX "" )
  #set_target_properties( pyPolyCSG PROPERTIES SUFFIX .so )
ELSEIF( UNIX )
  # Unix specific library naming
ELSEIF( WIN32 )
  # Windows specific library naming
ENDIF( APPLE )
set_target_properties( pyPolyCSG PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ..
  LIBRARY_OUTPUT_DIRECTORY_DEBUG ..
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ..
)

#add_executable( pyPolyCSG_test source/boolean_test.cpp )
#target_link_libraries( pyPolyCSG_test pyPolyCSG )
