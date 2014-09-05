// Copyright (C) 2012 Benjamin Kehlet
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2012-10-31
// Last changed: 2012-11-14

// Some utilities for working with cgal polyhedrons


#ifndef __POLYHEDRON_UTILS_H
#define __POLYHEDRON_UTILS_H

//#ifdef HAS_CGAL

#include <CGAL/basic.h>
#include<CGAL/Nef_polyhedron_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include<CGAL/Polyhedron_incremental_builder_3.h>
#include<CGAL/Polyhedron_3.h>
#include<CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Nef_polyhedron_3.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//#include <CGAL/Mesh_triangulation_3.h>
//#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
//#include <CGAL/Triangulation_vertex_base_with_info_3.h>
//#include <CGAL/Triangulation_cell_base_with_info_3.h>

//#include <CGAL/IO/Polyhedron_iostream.h>
//#include <CGAL/Bbox_3.h>

//#include <CGAL/Mesh_criteria_3.h>
//#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
//#include <CGAL/make_mesh_3.h>


//#include "cgal_csg3d.h"
//#include "self_intersect.h"

//// Criteria
//typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron;


namespace polyhcsg
{

//  class PolyhedronUtils
//  {
//  public:
//    static void readSurfaceFile(std::string filename,
//                                csg::Exact_Polyhedron_3& p);
//    static void readSTLFile(std::string filename, csg::Exact_Polyhedron_3& p);
//    static CGAL::Bbox_3 getBoundingBox(csg::Polyhedron_3& polyhedron);
//    static double getBoundingSphereRadius(csg::Polyhedron_3& polyhedron);

bool has_degenerate_facets(Polyhedron& p, double threshold);

void remove_degenerate_facets(Polyhedron& p, const double threshold);

//    template <typename Polyhedron>
//    bool has_self_intersections(Polyhedron& p)
//    {
//      typedef typename Polyhedron::Triangle_3 Triangle;
//      typedef typename std::back_insert_iterator<std::list<Triangle> > OutputIterator;
//
//      std::list<Triangle> triangles; // intersecting triangles
//      ::self_intersect<Polyhedron::Polyhedron_3, Polyhedron::Kernel,
//          OutputIterator>(p, std::back_inserter(triangles));
//
//      return triangles.size() > 0;
//    }
//  };
} // namespace polyhcsg

//#endif
#endif
