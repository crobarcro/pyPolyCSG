#include<vector>
#include<iterator>
#include<iostream>

#include"polyhedron.h"
#include"polyhedron_binary_op.h"

#if defined(CSG_USE_CGAL) && !defined(CSG_USE_CARVE)

#include "PolyhedronUtils.h"
//#include<CGAL/Nef_polyhedron_3.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include<CGAL/Polyhedron_incremental_builder_3.h>
//#include<CGAL/Polyhedron_3.h>
//#include<CGAL/IO/Polyhedron_iostream.h>
//#include <CGAL/Iterator_project.h>
//#include <CGAL/function_objects.h>

//#include <CGAL/Polyhedron_items_with_id_3.h>
//// Adaptor for Polyhedron_3
//#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
//// Simplification function
//#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
//// Visitor base
//#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
//// Stop-condition policy
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
//// Non-default cost and placement policies
//#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>



#elif defined(CSG_USE_CARVE)

#if defined(HAVE_CONFIG_H)
#include <carve/config.h>
#endif
#include <carve/interpolator.hpp>
#include <carve/csg_triangulator.hpp>
#include <carve/csg.hpp>

#endif

#if defined(CSG_USE_CGAL) && !defined(CSG_USE_CARVE)

//// the following is used for checking for degenerate meshes
//typedef Kernel::Point_3 Point ;
//typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Polyhedron_id;
//typedef Polyhedron_id::Halfedge_handle Halfedge_handle ;
//typedef Polyhedron_id::Vertex_handle Vertex_handle ;
//typedef CGAL::Surface_mesh_simplification::Edge_profile<Polyhedron_id> Profile ;
//typedef  boost::graph_traits<Polyhedron_id>::edges_size_type size_type ;
//
////*******************************************************************************************************************
//// -= stopping condition predicate =-
////
//// Determines whether the simplification has finished.
//// The arguments are (current_cost,vertex,vertex,is_edge,initial_pair_count,current_pair_count,surface) and the result is bool
////
////*******************************************************************************************************************
////
//// Stops when the cost is below some tolerance
////
//template<class ECM_>
//        class length_stop_predicate
//{
//        public:
//            typedef ECM_ ECM ;
//            typedef CGAL::Surface_mesh_simplification::Edge_profile<ECM> Profile ;
//            typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
//            typedef typename boost::graph_traits<ECM>::edges_size_type size_type ;
//            typedef typename CGAL::halfedge_graph_traits<ECM>::Point Point ;
//            typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel ;
//            typedef typename Kernel::FT FT ;
//        public :
//            length_stop_predicate( double tolerance ) : _tol(tolerance) {}
//            bool operator()( FT const& aCurrentCost // aCurrentCost
//                    , Profile const& //aEdgeProfile
//                    , size_type //aInitialCount
//                    , size_type //aCurrentCount
//                    ) const
//            {
//                std::cout << "sqrt aCurrentCost Is " << sqrt(CGAL::to_double(aCurrentCost)) << " _tol is " << _tol << std::endl;
//                return ( sqrt(CGAL::to_double(aCurrentCost)) > _tol );
//            }
//        private:
//            double _tol;
//};
//
// // The following is a Visitor that keeps track of the simplification process.
// // In this example the progress is printed real-time and a few statistics are
// // recorded (and printed in the end).
// //
// struct Stats
// {
//     Stats()
//     : collected(0)
//     , processed(0)
//     , collapsed(0)
//     , non_collapsable(0)
//     , cost_uncomputable(0)
//     , placement_uncomputable(0)
//     {}
//     std::size_t collected ;
//     std::size_t processed ;
//     std::size_t collapsed ;
//     std::size_t non_collapsable ;
//     std::size_t cost_uncomputable ;
//     std::size_t placement_uncomputable ;
// } ;
//
// struct My_visitor : CGAL::Surface_mesh_simplification::Edge_collapse_visitor_base<Polyhedron_id>
// {
//     My_visitor( Stats* s) : stats(s){}
//     // Called during the collecting phase for each edge collected.
//     void OnCollected( Profile const&, boost::optional<Kernel::FT> const& )
//     {
//         ++ stats->collected ;
//         std::cerr << "Edges collected: " << stats->collected << std::endl << std::flush ;
//     }
//
//     // Called during the processing phase for each edge selected.
//     // If cost is absent the edge won't be collapsed.
//     void OnSelected(Profile const&
//             ,boost::optional<Kernel::FT> cost
//             ,std::size_t initial
//             ,std::size_t current
//             )
//     {
//         ++ stats->processed ;
//         if ( !cost )
//             ++ stats->cost_uncomputable ;
//         if ( current == initial )
//             std::cerr << "\n" << std::flush ;
//         std::cerr << "\r" << current << std::flush ;
//     }
//
//     // Called during the processing phase for each edge being collapsed.
//     // If placement is absent the edge is left uncollapsed.
//     void OnCollapsing(Profile const&
//             ,boost::optional<Point> placement
//             )
//     {
//         if ( !placement )
//             ++ stats->placement_uncomputable ;
//     }
//
//     // Called for each edge which failed the so called link-condition,
//     // that is, which cannot be collapsed because doing so would
//     // turn the surface mesh into a non-manifold.
//     void OnNonCollapsable( Profile const& )
//     {
//         ++ stats->non_collapsable;
//     }
//
//     // Called AFTER each edge has been collapsed
//     void OnCollapsed( Profile const&, Vertex_handle )
//     {
//         ++ stats->collapsed;
//     }
//     Stats* stats ;
// } ;
#endif

namespace polyhcsg
{

#if defined(CSG_USE_CGAL) && !defined(CSG_USE_CARVE)

// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
    polyhedron t;

    polyhedron_builder( const polyhedron &p ){
        t = p;
    }
    void operator()( HDS& hds) {
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;


        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( t.num_vertices(), t.num_faces() );

        // add the polyhedron vertices
        for( int i=0; i<t.num_vertices(); i++ ){
            double x, y, z;
            t.get_vertex( i, x, y, z );
            B.add_vertex( Point( x, y, z ) );
        }

        std::vector<double> coords;
        std::vector<int> faces;
        t.output_store_in_mesh( coords, faces );

        int nverts, tmpi = 0;
        while( tmpi < (int)faces.size() ){
            nverts = faces[tmpi++];
            if( nverts != 3 )
                std::cout << "face has " << nverts << " vertices" << std::endl;
            B.begin_facet();
            for( int i=0; i<nverts; i++ ){
                B.add_vertex_to_facet( faces[tmpi++] );
            }
            B.end_facet();
        }

        // finish up the surface
        B.end_surface();
    }
};

Nef_polyhedron polyhedron_to_cgal( const polyhedron &p ){
    polyhedron tmp = p.triangulate();
    Polyhedron P;
    polyhedron_builder<HalfedgeDS> builder( tmp );
    P.delegate( builder );
    if( P.is_closed() )
        return Nef_polyhedron( P );
    else
        std::cout << "input polyhedron is not closed!" << std::endl;

    return Nef_polyhedron();
}


polyhedron cgal_to_polyhedron( const Nef_polyhedron &NP ){
//    Polyhedron_id P;
    Polyhedron P;
    polyhedron ret;

    if( NP.is_simple() )
    {
        NP.convert_to_polyhedron(P);

        // ensure there are no degenerate (too small) facets
        remove_degenerate_facets(P, 1.0e-8);

//        // ensure there are no degenerate (too small) edges
//
//        // The items in this polyhedron have an "id()" field
//        // which the default index maps used in the algorithm
//        // need to get the index of a vertex/edge.
//        // However, the Polyhedron_3 class doesn't assign any value to
//        // this id(), so we must do it here:
//        int index = 0 ;
//        for( Polyhedron_id::Halfedge_iterator eb = P.halfedges_begin()
//            , ee = P.halfedges_end()
//            ; eb != ee
//            ; ++ eb
//           )
//        {
//            eb->id() = index++;
//        }
//
//        index = 0 ;
//        for( Polyhedron_id::Vertex_iterator vb = P.vertices_begin()
//            , ve = P.vertices_end()
//            ; vb != ve
//            ; ++ vb
//           )
//        {
//            vb->id() = index++;
//        }
//
//        length_stop_predicate<Polyhedron_id> stop (0.001);
//        //CGAL::Surface_mesh_simplification::Count_ratio_stop_predicate<Polyhedron_id> stop(0.99);
//
//        Stats stats ;
//
//        My_visitor vis(&stats) ;
//
//        // The index maps are not explicitetly passed as in the previous
//        // example because the surface mesh items have a proper id() field.
//        // On the other hand, we pass here explicit cost and placement
//        // function which differ from the default policies, ommited in
//        // the previous example.
//        int r = CGAL::Surface_mesh_simplification::edge_collapse (
//                                    P
//                                    ,stop
//                                    ,CGAL::get_cost (CGAL::Surface_mesh_simplification::Edge_length_cost <Polyhedron_id>())
//                                    .get_placement(CGAL::Surface_mesh_simplification::Midpoint_placement<Polyhedron_id>())
//                                    .visitor (vis)
//                                                                 );
//
//        std::cout << "\nEdges collected: " << stats.collected
//                  << "\nEdges processed: " << stats.processed
//                  << "\nEdges collapsed: " << stats.collapsed
//                  << std::endl
//                  << "\nEdges not collapsed due to topological constraints: " << stats.non_collapsable
//                  << "\nEdge not collapsed due to cost computation constraints: " << stats.cost_uncomputable
//                  << "\nEdge not collapsed due to placement computation constraints: " << stats.placement_uncomputable
//                  << std::endl ;
//
//        std::cout << "\nFinished...\n" << r << " edges removed.\n"
//                  << (P.size_of_halfedges()/2) << " final edges.\n" ;


        std::vector<double> coords;
        std::vector<int> tris;
        int next_id = 0;
        std::map< Polyhedron::Vertex*, int > vid;
        for( Polyhedron::Vertex_iterator iter = P.vertices_begin(); iter != P.vertices_end(); iter++ )
        {
            coords.push_back( CGAL::to_double( (*iter).point().x() ) );
            coords.push_back( CGAL::to_double( (*iter).point().y() ) );
            coords.push_back( CGAL::to_double( (*iter).point().z() ) );
            vid[ &(*iter) ] = next_id++;
        }

        for( Polyhedron::Facet_iterator iter = P.facets_begin(); iter != P.facets_end(); iter++ )
        {
            Polyhedron::Halfedge_around_facet_circulator j = iter->facet_begin();
            tris.push_back( CGAL::circulator_size(j) );
            do
            {
                tris.push_back( std::distance(P.vertices_begin(), j->vertex()) );
            } while ( ++j != iter->facet_begin());
        }

        ret.initialize_load_from_mesh( coords, tris );
    }
    else
    {
        std::cout << "resulting polyhedron is not simple!" << std::endl;
    }
    return ret;
}

polyhedron polyhedron_union::operator()( const polyhedron &A, const polyhedron &B ){
    Nef_polyhedron a, b, c;
    try {
        a = polyhedron_to_cgal( A );
        b = polyhedron_to_cgal( B );
        c = (a + b).interior().closure();
        return cgal_to_polyhedron( c );
    } catch( std::exception &e ){
        return A;
    }
}

polyhedron polyhedron_difference::operator()( const polyhedron &A, const polyhedron &B ){
    Nef_polyhedron a, b, c;
    try {
        a = polyhedron_to_cgal( A );
        b = polyhedron_to_cgal( B );
        c = (a - b).interior().closure();
        return cgal_to_polyhedron( c );
    } catch( std::exception &e ){
        return A;
    }
}

polyhedron polyhedron_symmetric_difference::operator()( const polyhedron &A, const polyhedron &B ){
    Nef_polyhedron a, b, c;
    try {
        a = polyhedron_to_cgal( A );
        b = polyhedron_to_cgal( B );
        c = (a ^ b).interior().closure();
        return cgal_to_polyhedron( c );
    } catch( std::exception &e ){
        return A;
    }
}

polyhedron polyhedron_intersection::operator()( const polyhedron &A, const polyhedron &B ){
    Nef_polyhedron a, b, c;
    try {
        a = polyhedron_to_cgal( A );
        b = polyhedron_to_cgal( B );
        c = (a * b).interior().closure();
        return cgal_to_polyhedron( c );
    } catch( std::exception &e ){
        return A;
    }
}


#elif defined(CSG_USE_CARVE)

carve::mesh::MeshSet<3> *polyhedron_to_carve( const polyhedron &p ){
	std::vector<double> coords;
	std::vector<int> faces;

	p.output_store_in_mesh( coords, faces );

    std::vector<carve::mesh::MeshSet<3>::vertex_t*> v;
    std::vector<carve::mesh::MeshSet<3>::face_t *> f;
    for( int i=0; i<(int)coords.size(); i+=3 ){
        v.push_back( new carve::mesh::MeshSet<3>::vertex_t( carve::geom::VECTOR( coords[i+0], coords[i+1], coords[i+2]) ) );
    }
    int nverts, tmpi=0;
	while( tmpi < (int)faces.size() ){
		nverts = faces[tmpi++];
		std::vector<carve::mesh::MeshSet<3>::vertex_t*> face_verts;
		for( int i=0; i<nverts; i++ ){
			face_verts.push_back( v[faces[tmpi++]] );
		}
        carve::mesh::MeshSet<3>::face_t *tf = new carve::mesh::MeshSet<3>::face_t( face_verts.begin(), face_verts.end() );
		f.push_back( tf );
	}

	return new carve::mesh::MeshSet<3>( f );
}

polyhedron carve_to_polyhedron( carve::mesh::MeshSet<3> *p ){
	std::map< const carve::mesh::MeshSet<3>::vertex_t*, int > vid;
	std::vector<double> coords;
	std::vector<int> faces;

    int nextvid = 0;
    for( carve::mesh::MeshSet<3>::face_iter i=p->faceBegin(); i!=p->faceEnd(); ++i ){
        carve::mesh::MeshSet<3>::face_t *f = *i;

        std::vector<int> fvid;
        for (carve::mesh::MeshSet<3>::face_t::edge_iter_t e = f->begin(); e != f->end(); ++e) {
            carve::mesh::MeshSet<3>::vertex_t *tv = e->vert;
            if( vid.find(tv) == vid.end() ){
                vid[tv] = nextvid++;
                coords.push_back( tv->v.x );
                coords.push_back( tv->v.y );
                coords.push_back( tv->v.z );
            }
            fvid.push_back( vid[tv] );
        }
        faces.push_back( fvid.size() );
        for( int j=0; j<(int)fvid.size(); j++ ){
            faces.push_back( fvid[j] );
        }
    }

	polyhedron poly;
	poly.initialize_load_from_mesh( coords, faces );
	return poly;
}

polyhedron polyhedron_union::operator()( const polyhedron &A, const polyhedron &B ){
	carve::mesh::MeshSet<3> *pA = polyhedron_to_carve( A );
	carve::mesh::MeshSet<3> *pB = polyhedron_to_carve( B );
	carve::csg::CSG csg;
	carve::mesh::MeshSet<3> *pR = csg.compute( pA, pB, carve::csg::CSG::UNION );
	polyhedron R = carve_to_polyhedron( pR );
	delete pA;
	delete pB;
	delete pR;
	return R;
}

polyhedron polyhedron_difference::operator()( const polyhedron &A, const polyhedron &B ){
	carve::mesh::MeshSet<3> *pA = polyhedron_to_carve( A );
	carve::mesh::MeshSet<3> *pB = polyhedron_to_carve( B );
	carve::csg::CSG csg;
	carve::mesh::MeshSet<3> *pR = csg.compute( pA, pB, carve::csg::CSG::A_MINUS_B );
	polyhedron R = carve_to_polyhedron( pR );
	delete pA;
	delete pB;
	delete pR;
	return R;
}

polyhedron polyhedron_symmetric_difference::operator()( const polyhedron &A, const polyhedron &B ){
	carve::mesh::MeshSet<3> *pA = polyhedron_to_carve( A );
	carve::mesh::MeshSet<3> *pB = polyhedron_to_carve( B );
	carve::csg::CSG csg;
	carve::mesh::MeshSet<3> *pR = csg.compute( pA, pB, carve::csg::CSG::SYMMETRIC_DIFFERENCE );
	polyhedron R = carve_to_polyhedron( pR );
	delete pA;
	delete pB;
	delete pR;
	return R;
}

polyhedron polyhedron_intersection::operator()( const polyhedron &A, const polyhedron &B ){
	carve::mesh::MeshSet<3> *pA = polyhedron_to_carve( A );
	carve::mesh::MeshSet<3> *pB = polyhedron_to_carve( B );
	carve::csg::CSG csg;
	carve::mesh::MeshSet<3> *pR = csg.compute( pA, pB, carve::csg::CSG::INTERSECTION );
	polyhedron R = carve_to_polyhedron( pR );
	delete pA;
	delete pB;
	delete pR;
	return R;
}
#endif

bool poly_test ( polyhedron p)
{
    Nef_polyhedron NP = polyhedron_to_cgal( p );
    p = cgal_to_polyhedron( NP );
    return true;
}

} // namespace polyhcsg
