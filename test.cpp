#include "moab/Core.hpp"

using namespace moab;
using namespace std;

string test_file_name = "test.vtk";

int main( int argc, char** argv )
{
  Interface *iface = new Core;

    // need option handling here for input filename
 // if (argc > 1){
    //user has input a mesh file
 //   test_file_name = argv[1];
  //}  
    //load the mesh from vtk file
  ErrorCode rval = iface->load_mesh( "test_filename");
  //assert( rval == MB_SUCCESS);

    //get verts entities
  Range verts;
  rval = iface->get_entities_by_type(0, MBVERTEX, verts);
  //assert( rval == MB_SUCCESS);
    //get edge entities
  Range edges;
  rval = iface->get_entities_by_type(0, MBEDGE, edges);
  //assert(rval == MB_SUCCESS);

    //get triangular entities
  Range tri;
  rval = iface->get_entities_by_type(0, MBTRI, tri);
  //assert( rval == MB_SUCCESS);

    //get quad entities
  Range quads;
  rval = iface->get_entities_by_type(0, MBQUAD, quads);
  //assert(rval == MB_SUCCESS);

    //get hex entities
  Range hex;
  rval = iface->get_entities_by_type(0, MBHEX, hex);
  //assert(rval == MB_SUCCESS);

   //output the number of entities
  cout << "Number of vertices is " << verts.size() <<  endl;
  cout << "Number of edges is " << edges.size() <<  endl;
  cout << "Number of triangular faces is " << tri.size() <<  endl;
  cout << "Number of quad faces is " << quads.size() <<  endl;
  cout << "Number of hex is " << hex.size() <<  endl;
  
  return 0;
}