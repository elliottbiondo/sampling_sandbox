/*
 * Contents:  Random-number sampling using the Walker-Vose alias method,
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
 * M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
 * A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)
 */

#include "sampling.hpp"
#define MBI moab_instance()

Sampling *Sampling::instance_ = NULL;

void Sampling::create_instance(MBInterface *mb_impl)
{
  if (NULL == mb_impl) mb_impl = new MBCore();
  instance_ = new Sampling(mb_impl);
}

Sampling::Sampling(MBInterface *mb_impl)
   : mbImpl(mb_impl){}

Sampling::AliasTable::AliasTable(std::vector<double> p)
{
    n = p.size();
    prob.resize(n);
    alias.resize(n);
    std::vector<double> S(n);
    std::vector<double> L(n);
    double sum = 0;
    int i, a, g;

    // normalize
    for ( i=0; i<n; ++i ) {
        if( p[i]<0 ) {
            fprintf( stderr, "ransampl: invalid probability p[%i]<0\n", i );
        }
        sum += p[i];
    }
    if ( !sum ) {
        fprintf( stderr, "ransampl: no nonzero probability\n" );
    }
    for ( i=0; i<n; ++i )
        p[i] = (p[i] * n / sum);

    // Set separate index lists for small and large probabilities:
    int nS = 0;
    int nL = 0;
    for ( i=n-1; i>=0; --i ) {
        // at variance from Schwarz, we revert the index order
        if ( p[i]<1 )
            S[nS++] = i;
        else
            L[nL++] = i;
    }

    // Work through index lists
    while ( nS && nL ) {
        a = S[--nS]; // Schwarz's l
        g = L[--nL]; // Schwarz's g
        prob[a] = p[a];
        alias[a] = g;
        p[g] = p[g] + p[a] - 1;
        if ( p[g] < 1 )
            S[nS++] = g;
        else
            L[nL++] = g;
    }

    while ( nL )
        prob[ L[--nL] ] = 1;

    while ( nS )
        // can only happen through numeric instability
        prob[ S[--nS] ] = 1;
}

int Sampling::AliasTable::drawSample(double ran1, double ran2){
    int i = (int) n * ran1;
    return ran2 < prob[i] ? i : alias[i];
}


void Sampling::pdfFromMesh(char* fileName, char* tagName){

  MBEntityHandle loaded_file_set;
  // create meshset to load file into
  rval = MBI->create_meshset(MESHSET_SET, loaded_file_set );
  assert( rval == MB_SUCCESS );
  // load file
  rval = MBI->load_file( fileName, &loaded_file_set );
  assert( rval == MB_SUCCESS );
  // get entities
  rval = MBI->get_entities_by_dimension(loaded_file_set, 3, ves);
  int num_hex, num_quad;

  rval = MBI->get_number_entities_by_type(loaded_file_set, MBHEX, num_hex);
  rval = MBI->get_number_entities_by_type(loaded_file_set, MBQUAD, num_quad);
  if(num_hex = ves.size()) 
  {
    ve_type = MBHEX;
    verts_per_vol = 8;
  }
  else if (num_quad = ves.size())
  {
    ve_type = MBTET;
    verts_per_vol = 4;
  }
  else exit(1);

  std::cout << "results: "<< ve_type << std::endl;
 
  assert( rval == MB_SUCCESS );
  // get tag handle
  rval = MBI->tag_get_handle( tagName, moab::MB_TAG_VARLEN, MB_TYPE_DOUBLE, phtnSrcTag);
  // THIS ASSERT FAILS because we do not know number of energy groups a priori.
  //assert( rval == MB_SUCCESS );
  // get tag size
  int phtnSrcTagSize;
  rval = MBI->tag_get_bytes(phtnSrcTag, *(&phtnSrcTagSize));
  assert( rval == MB_SUCCESS );
  tagLen = phtnSrcTagSize/sizeof(double);

  phtnSrcData.resize(ves.size()*tagLen); 

  rval = MBI->tag_get_data( phtnSrcTag, ves, &phtnSrcData[0]);
  assert( rval == MB_SUCCESS );
  
  std::vector<double> volumes =  find_volumes();

  int i, j;
  for(i=0; i<ves.size(); ++i){
    for(j=0; j<tagLen; ++j){
     phtnSrcData[i*tagLen + j] *=  volumes[i];
    }
  }

}

std::vector<double> Sampling::find_volumes()
{
  std::vector<double> volumes (ves.size());
  int i;
  double volume;

   //std::vector<double> coords(3*pow(ves.size()+1,3));


    //MBRange verts;
    //rval = MBI()->get_entities_by_type(0, MBVERTEX, verts);
    //double x[geom_sets[0].size()];
    //double y[geom_sets[0].size()];
    //double z[geom_sets[0].size()];

    //rval = MBI()-> get_coords( verts, &x[0], &y[0], &z[0]);
    //if(gen::error(MB_SUCCESS!=result,"could not get coordinates of the vertices")) return result;

   //rval = MBI->get_entities_by_type(0, MBVERTEX, verts);
   //rval = MBI->get_coords(coords, coords.data());
   //rval = MBI->get_coords(connect.data(), connect.size(), coords.data());


  MBRange verts;
  rval = MBI->get_entities_by_type(0, MBVERTEX, verts);
  std::vector<double> x (verts.size());
  std::vector<double> y (verts.size());
  std::vector<double> z (verts.size());
  MBI->get_coords(verts, &x[0], &y[0], &z[0]);


//std::vector< std::vector< double > > coords ( verts.size(), std::vector<double> (3));
  std::vector< std::vector< double > > coords;

  for(i=0; i<verts.size(); i++){
    std::vector< double> row;
    row.push_back(x[i]);
    row.push_back(y[i]);
    row.push_back(z[i]);
    coords.push_back(row);
  }

   for(i=0; i<27; i++){
     std::cout << coords[i][0] << coords[i][1]<< coords[i][2]<<std::endl;
   }

  double a[3] = {4.4,5.5,6.6};
 //= {1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.1, 11.1, 12.2, 13.3, 14.4, 15.5, 16.6, 17.7, 18.8};
  volume = measure(ve_type, 8, a);
  std::cout << volume << std::endl;

  //for(i=0; i<ves.size(); ++i){
  //  volumes[i] = MBI->measure(ve_type, verts_per_vol, coords[i]);
 //` }
  return volumes;
  
}



int main(int argc, char* argv[])
{
  double my_array[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> my_vec(&my_array[0], &my_array[0]+5);

  Sampling::AliasTable myTable(my_vec);

  int answers[] = {0, 0, 0, 0, 0};
  int N = 1000000;
  int i = 0;
  double rand1, rand2;
  for(i=0; i<N; i++){
     rand1 = (double) rand()/RAND_MAX;
     rand2 = (double) rand()/RAND_MAX;
     answers[myTable.drawSample(rand1, rand2)]++;
  }

  printf("bin |  prob  | expected prob\n");
  for(i=0; i<5; i++){
    printf("%i    %f   %f \n", i+1, (double) answers[i]/N, (double) (i+1)/15.0);
  }

 Sampling& sampling = *Sampling::instance();
 sampling.pdfFromMesh(argv[1], argv[2]);

//  int j;
//  for(i=0; i<sampling.ves.size(); ++i){
//    for(j=0; j<sampling.tagLen; ++j){
//    std::cout << sampling.phtnSrcData[i*sampling.tagLen + j] << std::endl;
//    }
//  }

return 0;
}
