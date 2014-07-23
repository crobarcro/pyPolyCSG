#ifndef PYPOLYHEDRON_H
#define PYPOLYHEDRON_H

#include<iostream>
#include<stdexcept>
#include<boost/python.hpp>
#include "polyhedron.h"

/**
 @file   pypolyhedron.h
 @author James Gregson (james.gregson@gmail.com)
 @brief Wrapper for the polyhedron class which forms the basis for the CSG library.  Includes methods for exporting values in python formats.
*/

/**
 @brief pypolyhedron class, thin wrapper for polyhedron class, wth python specific methods.
*/
class pypolyhedron : public polyhcsg::polyhedron {

public:
	/**
	 @brief default constructor
	*/
	pypolyhedron();

	/**
	 @brief constructor from polyhedron object
	*/
	pypolyhedron( const polyhcsg::polyhedron &in );

	/**
	 @brief copy constructor, performs a deep copy of the data
	*/
	pypolyhedron( const pypolyhedron &in );


    /**
     @brief returns a tuple containing the vertex_id'th vertex's coordinates
     @param[in] vertex_id input vertex id
     @return tuple containing the vertex_id'th vertex coordinates
    */
    boost::python::tuple py_get_vertex_coordinates( int vertex_id );

    /**
     @brief returns a list containing the face_id'th vertex indices
     @param[in] face_id the id of the face to get the vertex indices of
     @return list of vertex indices for the face_id'th face
    */
    boost::python::list py_get_face_vertices( int face_id );

    /**
     @brief returns a numpy array of mesh vertex coordinates as a 2D numpy array
     @return numpy array of vertex coordinates
    */
    boost::python::numeric::array py_get_vertices();

    /**
     @brief temporarily triangulates the current pypolyhedron and fills a 2D numpy array with the triangle vertex indices
     @return numpy array of triangle vertex indices
    */
    boost::python::numeric::array py_get_triangles();

    /**
     @brief python version of mult_matrix_3 to get around boost::python argument count limits
    */
    pypolyhedron py_mult_matrix_3( const boost::python::list &m ) const;

    /**
     @brief python version of mult_matrix_4 to get around boost::pythoon argument count limits
    */
    pypolyhedron py_mult_matrix_4( const boost::python::list &m ) const;

};


/**
 @brief convenience function for loading a pypolyhedron from a mesh without first constructing an object
 @param[in] filename input mesh filename
 @return loaded pypolyhedron (if successful) or an empty pypolyhedron
*/
pypolyhedron load_mesh_file( const char *filename );

/**
 @brief convenience function for creating a sphere, basically wraps pypolyhedron::initialize_create_sphere()
 @param[in] radius sphere radius
 @param[in] is_centered set to true if the sphere should be centered at the origin, false otherwise
 @param[in] hsegments number of segments around the equator
 @param[in] vsegments number of segments pole-to-pole
 @return generated pypolyhedron if successful or an empty pypolyhedron otherwise
*/
pypolyhedron sphere( const double radius, const bool is_centered=false, const int hsegments=20, const int vsegments=20 );

/**
 @brief convenience function for creating a box, wraps pypolyhedron::initialize_create_box()
 @param[in] size_x size along the x-axis
 @param[in] size_y size along the y-axis
 @param[in] size_z size along the z-axis
 @param[in] is_centered set this to true if the box should be centered
 @return generated pypolyhedron if successful, or an empty pypolyhedron otherwise
*/
pypolyhedron box( const double size_x, const double size_y, const double size_z, const bool is_centered=false );

/**
 @brief convenience function for creating a cylinder, wraps pypolyhedron::initialize_create_cylinder()
 @param[in] radius cylinder radius
 @param[in] height cylinder height
 @param[in] is_centered set this to true if the cylinder should be centered on the origin
 @param[in] segments number of segments
 @return generated pypolyhedron if successful, or an empty pypolyhedron otherwise
*/
pypolyhedron cylinder( const double radius, const double height, const bool is_centered=false, const int segments=20 );

/**
 @brief convenience function for creating a cone, wraps pypolyhedron::initialize_create_cone()
 @param[in] radius radius of the base of the cone
 @param[in] height height of the cone
 @param[in] is_centered set this to true if the cone should be centered on the origin
 @param[in] segments number of segments to discretize the cone
 @return generated pypolyhedron if successful, or an empty pypolyhedron otherwise
*/
pypolyhedron cone( const double radius, const double height, const bool is_centered=false, const int segments=20 );

/**
 @brief convenience function for creating a torus, wraps pypolyhedron::initialize_create_torus()
 @param[in] radius_major the major (larger) radius of the torus
 @param[in] radius_minor the minor (smaller) radius of the torus
 @param[in] is_centered set this to true if the torus should be centered on the origin
 @param[in] major_segments number of segments to discretize the major radius
 @param[in] minor_segments number of segments to discretize the minor radius
 @return generated pypolyhedron or an empty pypolyhedron otherwise
*/
pypolyhedron torus( const double radius_major, const double radius_minor, const bool is_centered=false, const int major_segments=20, const int minor_segments=20 );

/**
 @brief generates an extrusion from the input 2D profile
 @param[in] coords coordinates of the profile vertices, as in initialize_create_extrude()
 @param[in] lines coordinate of the profile lines, as in initialize_create_extrude()
 @param[in] distance the distance to extrude the contour
 @return true if successful, false otherwise
 */
pypolyhedron extrusion( const std::vector<double> &coords, const std::vector<int> &lines, const double distance );

/**
 @brief generates a surface of revolution (SOR )from the input 2D profile by revolving a contour around the x-axis
 @param[in] coords coordinates of the profile vertices, as in initialize_create_extrude()
 @param[in] lines coordinate of the profile lines, as in initialize_create_extrude()
 @param[in] angle angle by which to rotate, to generate partial SORs
 @param[in] segments number of segments to create
 @return true if successful, false otherwise
 */
pypolyhedron surface_of_revolution( const std::vector<double> &coords, const std::vector<int> &lines, const double angle=360.0, const int segments=20 );

#endif
