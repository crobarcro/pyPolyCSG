cmake_minimum_required( VERSION 2.6 )

# create a basic project
project( polyhcsg )

set( BOOLEAN_SOURCES
  source/mesh_functions.cpp
  source/mesh_io.cpp
  source/polyhedron_binary_op.cpp
  source/polyhedron_unary_op.cpp
  source/polyhedron.cpp
  source/triangulate.cpp
)

set( BOOLEAN_HEADERS
  include/mesh_functions.h
  include/mesh_io.h
  include/polyhedron_binary_op.h
  include/polyhedron_unary_op.h
  include/polyhedron.h
  include/triangulate.h
)

IF( NOT CSG_NO_CARVE )
find_package( CARVE )
IF( CARVE_FOUND )
  message( STATUS "Found system Carve installation!" )
  set( INCLUDE_DIRS ${INCLUDE_DIRS} ${CARVE_INCLUDE_DIRS} )
  set( LIBS ${LIBS} ${CARVE_LIBRARIES} )
  add_definitions( -DCSG_USE_CARVE )
ENDIF( CARVE_FOUND )
ENDIF( NOT CSG_NO_CARVE )

if( NOT CSG_NO_CGAL )
  find_package( CGAL )
  IF( CGAL_FOUND )
    message( STATUS "Found system CGAL installation!" )
    set( INCLUDE_DIRS ${INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS} )
    set( LIBS ${LIBS} ${CGAL_LIBRARIES} )
    include( ${CGAL_USE_FILE} )
    add_definitions( -DCSG_USE_CGAL )
  ENDIF( CGAL_FOUND )
ENDIF( NOT CSG_NO_CGAL )

IF( NOT CARVE_FOUND AND NOT CGAL_FOUND )
  message( STATUS "No system Carve or CGAL installation found, falling back to Carve in third_party/" )
  SET( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_SOURCE_DIR}/third_party )
  find_package( CARVE )
  IF( CARVE_FOUND )
    message( STATUS "Found third_party/ Carve installation!" )
    set( INCLUDE_DIRS ${INCLUDE_DIRS} ${CARVE_INCLUDE_DIRS} )
    set( LIBS ${LIBS} ${CARVE_LIBRARIES} )
    add_definitions( -DCSG_USE_CARVE )
  ELSE(CARVE_FOUND)
    message( ERROR "One of CGAL or Carve must be installed to build pyPolyCSG!" )
  ENDIF(CARVE_FOUND)
ENDIF( NOT CARVE_FOUND AND NOT CGAL_FOUND )

include_directories( include ${INCLUDE_DIRS} )
add_library( polyhcsg SHARED ${BOOLEAN_SOURCES} ${BOOLEAN_HEADERS} )
target_link_libraries( polyhcsg ${LIBS} )
IF( APPLE )
  # OS-X specific library naming
  #set_target_properties( pyPolyCSG PROPERTIES PREFIX "" )
  #set_target_properties( pyPolyCSG PROPERTIES SUFFIX .so )
ELSEIF( UNIX )
  # Unix specific library naming
ELSEIF( WIN32 )
  # Windows specific library naming
ENDIF( APPLE )
set_target_properties( polyhcsg PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ..
  LIBRARY_OUTPUT_DIRECTORY_DEBUG ..
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ..
)
