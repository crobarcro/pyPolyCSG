#include<map>
#include<cmath>
#include<iostream>

#include"mesh_io.h"
#include"polyhedron.h"
#include"polyhedron_unary_op.h"
#include"polyhedron_binary_op.h"
#include"triangulate.h"

namespace polyhcsg {

polyhedron load_mesh_file( const char *filename ){
	polyhedron p;
	p.initialize_load_from_file( filename );
	return p;
}

polyhedron sphere( const double radius, const bool is_centered, const int hsegments, const int vsegments ){
	polyhedron p;
	p.initialize_create_sphere( radius, is_centered, hsegments, vsegments );
	return p;
}

polyhedron box( const double size_x, const double size_y, const double size_z, const bool is_centered ){
	polyhedron p;
	p.initialize_create_box( size_x, size_y, size_z, is_centered );
	return p;
}

polyhedron cylinder( const double radius, const double height, const bool is_centered, const int segments ){
	polyhedron p;
	p.initialize_create_cylinder( radius, height, is_centered, segments );
	return p;
}

polyhedron cone( const double radius, const double height, const bool is_centered, const int segments ){
	polyhedron p;
	p.initialize_create_cone( radius, height, is_centered, segments );
	return p;
}

polyhedron torus( const double radius_major, const double radius_minor, const bool is_centered, const int major_segments, const int minor_segments ){
	polyhedron p;
	p.initialize_create_torus( radius_major, radius_minor, is_centered, major_segments, minor_segments );
	return p;
}

polyhedron extrusion( const std::vector<double> &coords, const std::vector<int> &lines, const double distance ){
    polyhedron p;
    p.initialize_create_extrusion( coords, lines, distance );
    return p;
}

polyhedron surface_of_revolution( const std::vector<double> &coords, const std::vector<int> &lines, const double angle, const int segments ){
    polyhedron p;
    p.initialize_create_surface_of_revolution( coords, lines, angle, segments );
    return p;
}


// class methods

polyhedron::polyhedron(){
	m_coords.clear();
	m_faces.clear();
    m_faces_start.clear();
}

polyhedron::polyhedron( const polyhedron &in ){
	m_coords = in.m_coords;
	m_faces  = in.m_faces;
    m_faces_start = in.m_faces_start;
}

polyhedron::polyhedron( const polyhedron* in ){
	m_coords = in->getCoords ();
	m_faces  = in->getFaces ();
    m_faces_start = in->getFacesStart ();
}

std::vector<double> polyhedron::getCoords() const {
    return m_coords;
}

std::vector<int> polyhedron::getFaces() const {
    return m_faces;
}

std::vector<int> polyhedron::getFacesStart() const {
    return m_faces_start;
}

bool polyhedron::initialize_load_from_file( const char *filename ){
	// temporary polyhedron data
	std::vector<double> coords;
	std::vector<int>    faces;

	// attempt to load the mesh file
	if( !load_mesh_file( filename, coords, faces ) )
		return false;

	// mesh loaded successfully, now try and build a polyhedron
	return initialize_load_from_mesh( coords, faces );
}

bool polyhedron::initialize_load_from_mesh( const std::vector<double> &coords, const std::vector<int> &faces ){
	// TODO: add checks for self-intersection, non-manifold edges

	// for now, just copy the arrays over
	m_coords = coords;
	m_faces  = faces;

    // build the face_start array, which
    // gives the starting index of each face
    // (to the entry containing the number of
    // vertices in the face).
    int i=0;
    m_faces_start.clear ();
    while( i < m_faces.size() ){
        m_faces_start.push_back( i );
        i += m_faces[i]+1;
    }

	return true;
}

/*
 2013-03-03 - Fixed so number of vertical segments was correct
*/
bool polyhedron::initialize_create_sphere( const double radius, const bool is_centered, const int hsegments, const int vsegments ){
	std::vector<double> coords;
	std::vector<int>    faces;

	for( int j=1; j<vsegments; j++ ){
		double phi = M_PI*double(j)/double(vsegments);
		for( int i=0; i<hsegments; i++ ){
			double theta = 2.0*M_PI*double(i)/double(hsegments);
			double x = radius*sin(phi)*cos(-theta);
			double y = radius*sin(phi)*sin(-theta);
			double z = radius*cos(phi);
			if( !is_centered ){
				x += radius;
				y += radius;
				z += radius;
			}
			coords.push_back( x );
			coords.push_back( y );
			coords.push_back( z );
		}
	}
    int bottom_vtx = coords.size()/3;
	coords.push_back( is_centered ? 0.0 : radius );
	coords.push_back( is_centered ? 0.0 : radius );
	coords.push_back( is_centered ? radius : 2.0*radius );

    int top_vtx = coords.size()/3;
	coords.push_back( is_centered ? 0.0 : radius );
	coords.push_back( is_centered ? 0.0 : radius );
	coords.push_back( is_centered ? -radius : 0.0 );

	for( int j=0; j<vsegments-2; j++ ){
		for( int i=0; i<hsegments; i++ ){
			faces.push_back( 4 );
			faces.push_back( (i+0)%hsegments+((j+0)%vsegments)*hsegments );
			faces.push_back( (i+1)%hsegments+((j+0)%vsegments)*hsegments );
			faces.push_back( (i+1)%hsegments+((j+1)%vsegments)*hsegments );
			faces.push_back( (i+0)%hsegments+((j+1)%vsegments)*hsegments );
		}
	}

	for( int i=0; i<hsegments; i++ ){
		faces.push_back( 3 );
		faces.push_back( bottom_vtx );
		faces.push_back( (i+1)%hsegments );
		faces.push_back( i );
	}
	for( int i=0; i<hsegments; i++ ){
		faces.push_back( 3 );
		faces.push_back( top_vtx );
		faces.push_back( i+(vsegments-2)*hsegments );
		faces.push_back( (i+1)%hsegments+(vsegments-2)*hsegments );
	}

	return initialize_load_from_mesh( coords, faces );
}

bool polyhedron::initialize_create_box( const double size_x, const double size_y, const double size_z, const bool is_centered ){
	std::vector<double> coords;
	std::vector<int> faces;

	if( is_centered ){
		coords.push_back( -size_x/2. ); coords.push_back( -size_y/2. ); coords.push_back( -size_z/2. );
		coords.push_back(  size_x/2. ); coords.push_back( -size_y/2. ); coords.push_back( -size_z/2. );
		coords.push_back(  size_x/2. ); coords.push_back(  size_y/2. ); coords.push_back( -size_z/2. );
		coords.push_back( -size_x/2. ); coords.push_back(  size_y/2. ); coords.push_back( -size_z/2. );

		coords.push_back( -size_x/2. ); coords.push_back( -size_y/2. ); coords.push_back(  size_z/2. );
		coords.push_back(  size_x/2. ); coords.push_back( -size_y/2. ); coords.push_back(  size_z/2. );
		coords.push_back(  size_x/2. ); coords.push_back(  size_y/2. ); coords.push_back(  size_z/2. );
		coords.push_back( -size_x/2. ); coords.push_back(  size_y/2. ); coords.push_back(  size_z/2. );
	} else {
		coords.push_back( 0.0*size_x ); coords.push_back( 0.0*size_y ); coords.push_back( 0.0*size_z );
		coords.push_back( 1.0*size_x ); coords.push_back( 0.0*size_y ); coords.push_back( 0.0*size_z );
		coords.push_back( 1.0*size_x ); coords.push_back( 1.0*size_y ); coords.push_back( 0.0*size_z );
		coords.push_back( 0.0*size_x ); coords.push_back( 1.0*size_y ); coords.push_back( 0.0*size_z );

		coords.push_back( 0.0*size_x ); coords.push_back( 0.0*size_y ); coords.push_back( 1.0*size_z );
		coords.push_back( 1.0*size_x ); coords.push_back( 0.0*size_y ); coords.push_back( 1.0*size_z );
		coords.push_back( 1.0*size_x ); coords.push_back( 1.0*size_y ); coords.push_back( 1.0*size_z );
		coords.push_back( 0.0*size_x ); coords.push_back( 1.0*size_y ); coords.push_back( 1.0*size_z );
	}

	faces.push_back( 4 );
	faces.push_back( 0 );
	faces.push_back( 1 );
	faces.push_back( 5 );
	faces.push_back( 4 );

	faces.push_back( 4 );
	faces.push_back( 2 );
	faces.push_back( 3 );
	faces.push_back( 7 );
	faces.push_back( 6 );

	faces.push_back( 4 );
	faces.push_back( 3 );
	faces.push_back( 2 );
	faces.push_back( 1 );
	faces.push_back( 0 );

	faces.push_back( 4 );
	faces.push_back( 4 );
	faces.push_back( 5 );
	faces.push_back( 6 );
	faces.push_back( 7 );

	faces.push_back( 4 );
	faces.push_back( 1 );
	faces.push_back( 2 );
	faces.push_back( 6 );
	faces.push_back( 5 );

	faces.push_back( 4 );
	faces.push_back( 3 );
	faces.push_back( 0 );
	faces.push_back( 4 );
	faces.push_back( 7 );

	return initialize_load_from_mesh( coords, faces );
}

bool polyhedron::initialize_create_cylinder( const double radius, const double height, const bool is_centered, const int segments ){
	std::vector<double> coords;
	std::vector<int>    faces;

	double dtheta = M_PI*2.0/double(segments);
	if( is_centered ){
		for( int i=0; i<segments; i++ ){
			double c = cos(-i*dtheta);
			double s = sin(-i*dtheta);
			coords.push_back( radius*c );
			coords.push_back( -height/2.0 );
			coords.push_back( radius*s );
		}
	} else {
		for( int i=0; i<segments; i++ ){
			double c = cos(-i*dtheta);
			double s = sin(-i*dtheta);
			coords.push_back( radius*c+radius );
			coords.push_back( 0.0 );
			coords.push_back( radius*s+radius );
		}
	}
	for( int i=0; i<segments*3; i+=3 ){
		coords.push_back( coords[i+0] );
		coords.push_back( coords[i+1]+height );
		coords.push_back( coords[i+2] );
	}

	faces.push_back( segments );
	for( int i=0; i<segments; i++ ){
		faces.push_back( segments-i-1 );
	}
	faces.push_back( segments );
	for( int i=0; i<segments; i++ ){
		faces.push_back( i+segments );
	}
	for( int i=0; i<segments; i++ ){
		faces.push_back( 4 );
		faces.push_back( i );
		faces.push_back( (i+1)%segments );
		faces.push_back( (i+1)%segments+segments );
		faces.push_back( i+segments );
	}

	return initialize_load_from_mesh( coords, faces );
}


bool polyhedron::initialize_create_cone( const double radius, const double height, const bool is_centered, const int segments ){
	std::vector<double> coords;
	std::vector<int>    faces;

	for( int i=0; i<segments; i++ ){
		double theta = -2.0*M_PI*double(i)/double(segments);
		if( is_centered ){
			coords.push_back( radius*cos(theta) );
			coords.push_back( -height/2.0 );
			coords.push_back( radius*sin(theta) );
		} else {
			coords.push_back( radius*cos(theta)+radius );
			coords.push_back( 0.0 );
			coords.push_back( radius*sin(theta)+radius );
		}
	}
	if( is_centered ){
		coords.push_back( 0.0 );
		coords.push_back( height/2.0 );
		coords.push_back( 0.0 );
	} else {
		coords.push_back( radius );
		coords.push_back( height );
		coords.push_back( radius );
	}

	for( int i=0; i<segments; i++ ){
		faces.push_back( 3 );
		faces.push_back( i );
		faces.push_back( (i+1)%segments );
		faces.push_back( segments );
	}
	faces.push_back( segments );
	for( int i=segments-1; i>=0; i-- ){
		faces.push_back( i );
	}

	return initialize_load_from_mesh( coords, faces );
}

bool polyhedron::initialize_create_torus( const double radius_major, const double radius_minor, const bool is_centered, const int major_segments, const int minor_segments ){
	std::vector<double> coords;
	std::vector<int>    faces;

	for( int j=0; j<major_segments; j++ ){
		double phi = 2.0*M_PI*double(j)/double(major_segments);
		for( int i=0; i<minor_segments; i++ ){
			double theta = 2.0*M_PI*double(i)/double(minor_segments);
			if( is_centered ){
				coords.push_back( (radius_major+radius_minor*cos(theta))*cos(phi) );
				coords.push_back( radius_minor*sin(theta) );
				coords.push_back( (radius_major+radius_minor*cos(theta))*sin(phi) );
			} else {
				coords.push_back( (radius_major+radius_minor*cos(theta))*cos(phi)+radius_minor+radius_major );
				coords.push_back( radius_minor*sin(theta)+radius_minor );
				coords.push_back( (radius_major+radius_minor*cos(theta))*sin(phi)+radius_minor+radius_major );
			}
		}
	}
	for( int j=0; j<major_segments; j++ ){
		for( int i=0; i<minor_segments; i++ ){
			faces.push_back( 4 );
			faces.push_back( (i+0)%minor_segments+((j+0)%major_segments)*minor_segments );
			faces.push_back( (i+1)%minor_segments+((j+0)%major_segments)*minor_segments );
			faces.push_back( (i+1)%minor_segments+((j+1)%major_segments)*minor_segments );
			faces.push_back( (i+0)%minor_segments+((j+1)%major_segments)*minor_segments );
		}
	}

	return initialize_load_from_mesh( coords, faces );
}

bool polyhedron::initialize_create_extrusion( const std::vector<double> &coords, const std::vector<int> &lines, const double distance ){
	std::vector<double> tcoords;
	std::vector<int>    tfaces;
	for( int i=0; i<(int)coords.size(); i+=2 ){
		tcoords.push_back( coords[i+0] );
		tcoords.push_back( coords[i+1] );
		tcoords.push_back( 0.0 );
	}
	for( int i=0; i<(int)coords.size(); i+=2 ){
		tcoords.push_back( coords[i+0] );
		tcoords.push_back( coords[i+1] );
		tcoords.push_back( distance );
	}

	tfaces.push_back( lines.size() );
	for( int i=0; i<(int)lines.size(); i++ ){
		tfaces.push_back( lines[lines.size()-i-1] );
	}
	tfaces.push_back( lines.size() );
	for( int i=0; i<(int)lines.size(); i++ ){
		tfaces.push_back( lines[i]+coords.size()/2 );
	}

	for( int i=0; i<(int)lines.size(); i++ ){
		tfaces.push_back( 4 );
		tfaces.push_back( lines[i] );
		tfaces.push_back( lines[(i+1)%lines.size()] );
		tfaces.push_back( lines[(i+1)%lines.size()]+coords.size()/2 );
		tfaces.push_back( lines[i]+coords.size()/2 );
	}

	return initialize_load_from_mesh( tcoords, tfaces );
}

bool polyhedron::initialize_create_extrusion( const std::vector<double> &coords, const std::vector<int> &lines, const double distance, const int segments, const double dtheta )
{
	typedef std::pair<int,int> ii_pair;

	std::vector<double> tcoords;
	std::vector<int>    tfaces;

	if( fabs(dtheta) > 45 ){
	    // large angles will cause problems
	}

	int next_id = 0;

    // convert rotation angle to radians
	double raddtheta = dtheta * M_PI / 180.0;
    double dZ = distance / segments;

	// ****          Here we build up the coordinates         ***** //

	// map of pairs of vertex ids making up the edges of the polyhedron
	// there will be one edge for every segment and line, i.e. a total
	// of segments * lines.size() edges
	std::map< ii_pair, int > vid_map;
	for( int segn=0; segn <= segments; segn++ )
	{
        // loop through each incremental rotation
		double c = cos(segn*raddtheta);
		double s = sin(segn*raddtheta);

		for( int linen=0; linen<(int)lines.size(); linen++ )
		{
            double x = coords[lines[linen]*2+0];
            double y = coords[lines[linen]*2+1];

		    // loop through each line rotating the coordinates
            tcoords.push_back( x*c - y*s ); // x-coord rotated
            tcoords.push_back( y*c + x*s ); // y-coord rotated
            tcoords.push_back( segn * dZ ); // z-coord calculated
            // add the edge to the edge list
            vid_map[ ii_pair(segn,linen) ] = next_id++;
 		}
	}

    // ****            Here we generate the faces            ***** //
	for( int segn=0; segn<segments; segn++ )
	{
	    // loop through each incremental rotation
		for( int linen=0; linen<(int)lines.size(); linen++ )
		{
		    // get index of the first vertex of the line from this segment
			int v0 = vid_map[ ii_pair((segn+0) % (segments+1), linen+0) ];
			// get index of the first vertex of the line from the next segment
			int v1 = vid_map[ ii_pair((segn+1) % (segments+1), linen+0) ];
			// get index of the second vertex of the same line from this segment
			int v2 = vid_map[ ii_pair((segn+0) % (segments+1), (linen+1) % lines.size()) ];
			// get index of the second vertex of the line from the next segment
			int v3 = vid_map[ ii_pair((segn+1) % (segments+1), (linen+1) % lines.size()) ];

            // create the quad face
            tfaces.push_back( 4 );
            tfaces.push_back( v1 );
            tfaces.push_back( v3 );
            tfaces.push_back( v2 );
            tfaces.push_back( v0 );
		}
	}

	// add caps for the top and bottom
    tfaces.push_back( lines.size() );
    for( int linen=0; linen<(int)lines.size(); linen++ )
    {
        tfaces.push_back( vid_map[ ii_pair(0,linen) ] );
    }
    tfaces.push_back( lines.size() );
    for( int linen=0; linen<(int)lines.size(); linen++ )
    {
        tfaces.push_back( vid_map[ ii_pair(segments,lines.size()-linen-1) ] );
    }

	return initialize_load_from_mesh( tcoords, tfaces );
}

bool polyhedron::initialize_create_surface_of_revolution( const std::vector<double> &coords, const std::vector<int> &lines, const double angle, const int segments ){
	typedef std::pair<int,int> ii_pair;

	std::vector<double> tcoords;
	std::vector<int>    tfaces;

	int lastseg = segments+1;
	if( fabs(angle-360.0) < 1e-5 ){
	    // below this tolerance the surface will be a closed loop
	    // so the last segment will actually be the first segment
		lastseg = segments;
	}

	int next_id = 0;

    // rotation angle between each rotated section
	double dtheta = angle*M_PI / (double(segments)*180.0);

	// ****          Here we build up the coordinates         ***** //

	// map of pairs of vertex ids making up the edges of the polyhedron
	// there will be one edge for every segment and line, i.e. a total
	// of segments * lines.size() edges
	std::map< ii_pair, int > vid_map;
	for( int segn=0; segn < lastseg; segn++ )
	{
        // loop through each incremental rotation
		double c = cos(segn*dtheta);
		double s = sin(segn*dtheta);

		for( int linen=0; linen<(int)lines.size(); linen++ )
		{
		    // loop through each line making up the shape to be rotated

		    // Test if the y coordinate of the second point making up the
		    // current line is zero to within a tolerance
			if( fabs(coords[ lines[linen]*2+1 ]) < 1e-6 )
			{
				if( segn == 0 )
				{
                    // This is the first segment so we allow the coordinate to
                    // be added, for other segments we will link back to this
                    // coordinate
					tcoords.push_back( coords[lines[linen]*2+0] ); // x-coord stays the same
					tcoords.push_back( c * coords[lines[linen]*2+1] ); // y-coord modified
					tcoords.push_back( s * coords[lines[linen]*2+1] ); // z-coord calculated
					// add the edge to the edge list
					vid_map[ ii_pair(segn,linen) ] = next_id++;
				}
				else
				{
                    // this is not the first segment, so instead we point to
                    // the location where this coordinate was first added, and
                    // don't add a new coordinate
					vid_map[ ii_pair(segn,linen) ] = vid_map[ ii_pair(segn-1,linen) ];
				}
			}
			else
			{
                // normal coordinate, not lying on y-axis
				tcoords.push_back( coords[lines[linen]*2+0] ); // x-coord stays the same
				tcoords.push_back( c*coords[lines[linen]*2+1] ); // y-coord modified
				tcoords.push_back( s*coords[lines[linen]*2+1] ); // z-coord calculated
				// add the edge to the edge list
				vid_map[ ii_pair(segn,linen) ] = next_id++;
			}
 		}
	}

    // If the angle is a full revolution, we must point to the first segment in
    // the vid map so the last set of faces is created, completing the surface
	if( lastseg == segments )
	{
        for( int linen=0; linen<(int)lines.size(); linen++ )
		{
            vid_map[ ii_pair(lastseg,linen) ] = vid_map[ ii_pair(0,linen) ];
        }
	}

    // ****            Here we generate the faces            ***** //
	for( int segn=0; segn<segments; segn++ )
	{
	    // loop through each incremental rotation
		for( int linen=0; linen<(int)lines.size(); linen++ )
		{
		    // get index of the first vertex of the line from this segment
			int v0 = vid_map[ ii_pair((segn+0) % lastseg, linen+0) ];
			// get index of the first vertex of the line from the next segment
			int v1 = vid_map[ ii_pair((segn+1) % lastseg, linen+0) ];
			// get index of the second vertex of the same line from this segment
			int v2 = vid_map[ ii_pair((segn+0) % lastseg, (linen+1) % lines.size()) ];
			// get index of the second vertex of the line from the next segment
			int v3 = vid_map[ ii_pair((segn+1) % lastseg, (linen+1) % lines.size()) ];

			if( v0 == v1 && v2 == v3 )
			{
				// do nothing, line-segment lies along axis, region contains no volume
			}
			else if( v0 == v1 )
			{
				// first point on axis, second off axis, triangular face
				tfaces.push_back( 3 );
				tfaces.push_back( v3 );
				tfaces.push_back( v2 );
				tfaces.push_back( v0 );
			}
			else if( v2 == v3 )
			{
				// first point off axis, second point on, triangular face
				tfaces.push_back( 3 );
				tfaces.push_back( v0 );
				tfaces.push_back( v1 );
				tfaces.push_back( v2 );
			}
			else
			{
				// both points off axis, quad face
				tfaces.push_back( 4 );
				tfaces.push_back( v1 );
				tfaces.push_back( v3 );
				tfaces.push_back( v2 );
				tfaces.push_back( v0 );
			}
		}
	}

	// if the angle is not a full revolution, add caps for the exposed sections
	if( lastseg != segments )
	{
		tfaces.push_back( lines.size() );
		for( int linen=0; linen<(int)lines.size(); linen++ )
		{
			tfaces.push_back( vid_map[ ii_pair(0,linen) ] );
		}
		tfaces.push_back( lines.size() );
		for( int linen=0; linen<(int)lines.size(); linen++ )
		{
			tfaces.push_back( vid_map[ ii_pair(segments,lines.size()-linen-1) ] );
		}
	}

	return initialize_load_from_mesh( tcoords, tfaces );
}

bool polyhedron::output_store_in_mesh( std::vector<double> &coords, std::vector<int> &faces ) const {
	coords = m_coords;
	faces  = m_faces;
	return true;
}

bool polyhedron::output_store_in_file( const char *filename ) const {
	return save_mesh_file( m_coords, m_faces, filename );
}

polyhedron polyhedron::triangulate() const {
	std::vector<double> coords;
	std::vector<int> faces;

	// the coordinate array will be the same
	coords = m_coords;

	int tmpi = 0;
	while( tmpi < (int)m_faces.size() ){
		int nverts = m_faces[tmpi];
		// if there are three vertices, just add them to the output
		if( nverts == 3 ){
			faces.push_back( m_faces[tmpi] );
			faces.push_back( m_faces[tmpi+1] );
			faces.push_back( m_faces[tmpi+2] );
			faces.push_back( m_faces[tmpi+3] );
		} else {
			// otherwise triangulate the face
			std::vector<int> tfaces;
			bool res = triangulate_simple_polygon( coords, &m_faces[tmpi], tfaces );
			if( !res ){
				std::cout << "failed to triangulate polygon with " << nverts << " vertices" << std::endl;
				faces.push_back( m_faces[tmpi] );
				for( int i=0; i<nverts; i++ ){
					faces.push_back( m_faces[tmpi+1+i] );
				}
			} else {
				for( int i=0; i<(int)tfaces.size(); i++ ){
					faces.push_back( tfaces[i] );
				}
			}
		}

		tmpi += nverts+1;
	}

	polyhedron poly;
	poly.initialize_load_from_mesh( coords, faces );
	return poly;
}

int polyhedron::num_vertices(){
    return m_coords.size()/3;
}

int polyhedron::num_faces(){
    return m_faces_start.size();
}

int polyhedron::num_face_vertices( int face_id ){
    if( face_id < 0 || face_id >= num_faces() ){
        throw std::range_error("invalid face id");
    }
    return m_faces[m_faces_start[face_id]];
}

void polyhedron::get_face_vertices( int face_id, int *vertex_id_list ){
    if( face_id < 0 || face_id >= num_faces() ){
        throw std::range_error("invalid face id");
    }
    int start = m_faces_start[face_id]+1;
    int n = m_faces[ m_faces_start[face_id] ];
    for( int i=0; i<n; i++ ){
        vertex_id_list[i] = m_faces[start+i];
    }
}

polyhedron polyhedron::translate( const double x, const double y, const double z ) const {
	polyhedron_translate op( x, y, z );
	return op( *this );
}

polyhedron polyhedron::rotate( const double theta_x, const double theta_y, const double theta_z ) const {
	polyhedron_rotate op( theta_x, theta_y, theta_z );
	return op( *this );
}

polyhedron polyhedron::scale( const double x, const double y, const double z ) const {
	polyhedron_scale op( x, y, z );
	return op( *this );
}

polyhedron polyhedron::mult_matrix_4(
                        double xx, double xy, double xz, double xa,
                        double yx, double yy, double yz, double ya,
                        double zx, double zy, double zz, double za,
                        double ax, double ay, double az, double aa ) const {
    polyhedron_multmatrix4 op( xx, xy, xz, xa,
                               yx, yy, yz, ya,
                               zx, zy, zz, za,
                               ax, ay, az, aa );

    return op( *this );
}


polyhedron polyhedron::mult_matrix_3(
        double xx,double xy,double xz,
        double yx,double yy,double yz,
        double zx,double zy,double zz ) const {
    polyhedron_multmatrix3 op( xx, xy, xz,
                               yx, yy, yz,
                               zx, zy, zz );

    return op( *this );
}


polyhedron polyhedron::operator+( const polyhedron &in ) const {
	polyhedron_union op;
	return op( *this, in );
}

polyhedron polyhedron::operator-( const polyhedron &in ) const {
	polyhedron_difference op;
	return op( *this, in );
}

polyhedron polyhedron::operator^( const polyhedron &in ) const {
	polyhedron_symmetric_difference op;
	return op( *this, in );
}

polyhedron polyhedron::operator*( const polyhedron &in ) const {
	polyhedron_intersection op;
	return op( *this, in );
}

} // namespace polyhcsg
