#include <vector>
#include "polyhcsg/polyhedron.h"

using namespace polyhcsg;

int main( int argc, const char* argv[] )
{
    // test extrude rotate
    std::vector<double> coords = {0.0, 0.0, 1.0, 0.0, 0.5, 1.0 };
    std::vector<int> lines = {0,1,2};
    double distance = 1;
    int segments = 2;
    double dTheta = 10;

    polyhedron p5;

    p5.initialize_create_extrusion( coords, lines, distance, segments, dTheta );

}
