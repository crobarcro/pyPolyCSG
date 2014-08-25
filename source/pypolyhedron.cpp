#include<map>
#include<cmath>
#include<iostream>

#include"pypolyhedron.h"


// convenience functions

pypolyhedron load_mesh_file( const char *filename ){
	pypolyhedron p;
	p.initialize_load_from_file( filename );
	return p;
}

pypolyhedron sphere( const double radius, const bool is_centered, const int hsegments, const int vsegments ){
	pypolyhedron p;
	p.initialize_create_sphere( radius, is_centered, hsegments, vsegments );
	return p;
}

pypolyhedron box( const double size_x, const double size_y, const double size_z, const bool is_centered ){
	pypolyhedron p;
	p.initialize_create_box( size_x, size_y, size_z, is_centered );
	return p;
}

pypolyhedron cylinder( const double radius, const double height, const bool is_centered, const int segments ){
	pypolyhedron p;
	p.initialize_create_cylinder( radius, height, is_centered, segments );
	return p;
}

pypolyhedron cone( const double radius, const double height, const bool is_centered, const int segments ){
	pypolyhedron p;
	p.initialize_create_cone( radius, height, is_centered, segments );
	return p;
}

pypolyhedron torus( const double radius_major, const double radius_minor, const bool is_centered, const int major_segments, const int minor_segments ){
	pypolyhedron p;
	p.initialize_create_torus( radius_major, radius_minor, is_centered, major_segments, minor_segments );
	return p;
}

pypolyhedron extrusion( const std::vector<double> &coords, const std::vector<int> &lines, const double distance ){
    pypolyhedron p;
    p.initialize_create_extrusion( coords, lines, distance );
    return p;
}

pypolyhedron surface_of_revolution( const std::vector<double> &coords, const std::vector<int> &lines, const double angle, const int segments ){
    pypolyhedron p;
    p.initialize_create_surface_of_revolution( coords, lines, angle, segments );
    return p;
}

// class methods

pypolyhedron::pypolyhedron(){

}

pypolyhedron::pypolyhedron( const polyhcsg::polyhedron &in ) : polyhedron() {
	m_coords = in.getCoords ();
	m_faces  = in.getFaces ();
    m_faces_start = in.getFacesStart();
}

pypolyhedron::pypolyhedron( const pypolyhedron &in ){
	m_coords = in.m_coords;
	m_faces  = in.m_faces;
    m_faces_start = in.m_faces_start;
}

boost::python::tuple pypolyhedron::py_get_vertex_coordinates( int vertex_id ){
    if( vertex_id < 0 || vertex_id >= num_vertices() ){
        throw std::range_error("invalid vertex id");
    }
    vertex_id *= 3;
    return boost::python::make_tuple( m_coords[vertex_id], m_coords[vertex_id+1], m_coords[vertex_id+2] );
}

boost::python::list pypolyhedron::py_get_face_vertices( int face_id ){
    if( face_id < 0 || face_id >= num_faces() ){
        throw std::range_error("invalid face id");
    }
    int start = m_faces_start[face_id]+1;
    int n = m_faces[ m_faces_start[face_id]];
    boost::python::list ret;
    for( int i=0; i<n; i++ ){
        ret.append( m_faces[start+i] );
    }
    return ret;
}

boost::python::numeric::array pypolyhedron::py_get_vertices(){
    boost::python::list tmp;
    for( int i=0; i<num_vertices(); i++ ){
        tmp.append( boost::python::make_tuple( m_coords[i*3+0], m_coords[i*3+1], m_coords[i*3+2] ) );
    }
    return boost::python::numeric::array( tmp );
}

boost::python::numeric::array pypolyhedron::py_get_triangles(){
    polyhedron tri = triangulate();
    boost::python::list tmp;
    for( int i=0; i<tri.num_faces(); i++ ){
        int vtx_id[3];
        tri.get_face_vertices( i, vtx_id );
        tmp.append( boost::python::make_tuple( vtx_id[0], vtx_id[1], vtx_id[2] ) );
    }
    return boost::python::numeric::array( tmp );
}

pypolyhedron pypolyhedron::py_mult_matrix_3( const boost::python::list &m ) const {
    double a[9];
    if( boost::python::len(m) != 9 )
        throw "expected 9 matrix coefficients";
    for( int i=0; i<9; i++ ){
        a[i] = boost::python::extract<double>( m[i] );
    }
    return pypolyhedron (mult_matrix_3( a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8] ));
}

pypolyhedron pypolyhedron::py_mult_matrix_4( const boost::python::list &m ) const {
    double a[16];
    if( boost::python::len(m) != 16 )
        throw "expected 16 matrix coefficients";
    for( int i=0; i<16; i++ ){
        a[i] = boost::python::extract<double>( m[i] );
    }
    return pypolyhedron (mult_matrix_4( a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15] ));
}
