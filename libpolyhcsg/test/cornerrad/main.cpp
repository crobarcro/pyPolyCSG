
#include "polyhcsg/polyhedron.h"
#include "polyhcsg/polyhedron_binary_op.h"

#include <stdio.h>

using namespace polyhcsg;

int main( int argc, const char* argv[] )
{
//	p = csg.polyhedron ();
//
//  p.makebox (0.6,0.6,0.6,0);

    polyhedron p1 = box( 0.6, 0.6, 0.6, false );

//  p2 = csg.polyhedron ();
//
//  p2.makesphere (0.5, 1, 37);

    polyhedron p2 = sphere( 0.5, true, 37 );

    polyhedron p_cornerrad;

    p_cornerrad = p1 - p2;

    polyhedron p3 = box( 0.6, 0.6, 0.6, false );

    polyhedron p4 = p3 - p_cornerrad;

}
