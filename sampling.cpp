/*
 * Contents:  Random-number sampling using the Walker-Vose alias method,
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
 * M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
 * A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)
 */


#include "sampling.hpp"
#define MBI moab_instance()

Sampling *Sampling::instance_ = NULL;
//MBInterface *MBI;

Sampling::AliasTable::AliasTable(std::vector<double> p){
    
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

void Sampling::create_instance(MBInterface *mb_impl)
{
  if (NULL == mb_impl) mb_impl = new MBCore();
  instance_ = new Sampling(mb_impl);
}

Sampling::Sampling(MBInterface *mb_impl)
   : mbImpl(mb_impl)
{
}


void Sampling::blash(){

 std::string input_filename="test.h5m";

 MBEntityHandle loaded_file_set;
 MBErrorCode rval;

 // create meshset to load file into
 rval = MBI->create_meshset(MESHSET_SET, loaded_file_set );
 assert( rval == MB_SUCCESS );
 // load file
 rval = MBI->load_file( input_filename.c_str(), &loaded_file_set );
 assert( rval == MB_SUCCESS );
 std::cout << "file loaded" << std::endl;

  MBRange hex;
  rval = MBI->get_entities_by_type( 0, MBHEX, hex );
  MBTag idxTag;
  MBTag phtnSrcTag;
  MBTag phtnSrcTag2;

  rval = MBI->tag_get_handle( "idx", moab::MB_TAG_VARLEN, MB_TYPE_INTEGER, idxTag);
  rval = MBI->tag_get_handle( "phtn_src", moab::MB_TAG_VARLEN, MB_TYPE_DOUBLE, phtnSrcTag);
  rval = MBI->tag_get_handle( "phtn_src2", moab::MB_TAG_VARLEN, MB_TYPE_DOUBLE, phtnSrcTag2);

  int idxTagSize;
  int phtnSrcTagSize;
  int phtnSrcTagSize2;

  rval = MBI->tag_get_bytes(idxTag, *(&idxTagSize));
  rval = MBI->tag_get_bytes(phtnSrcTag, *(&phtnSrcTagSize));
  rval = MBI->tag_get_bytes(phtnSrcTag2, *(&phtnSrcTagSize2));

  std::vector<int> idxData;
  idxData.resize(hex.size()*idxTagSize/sizeof(int)); 

  std::vector<double> phtnSrcData;
  phtnSrcData.resize(hex.size()*phtnSrcTagSize/sizeof(double)); 

  std::vector<double> phtnSrcData2;
  phtnSrcData2.resize(hex.size()*phtnSrcTagSize2/sizeof(double)); 

  rval = MBI->tag_get_data( idxTag, hex, &idxData[0]);
  rval = MBI->tag_get_data( phtnSrcTag, hex, &phtnSrcData[0]);
  rval = MBI->tag_get_data( phtnSrcTag2, hex, &phtnSrcData2[0]);
  
  int i;
  for(i=0; i<hex.size(); ++i){
    std::cout << idxData[i] <<" "<< phtnSrcData[i] <<" "<<phtnSrcData2[2*i]<<" "<< phtnSrcData2[2*i+1] <<" "<< std::endl;
  }
  
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
 sampling.blash();
return 0;
}

//MBInterface *MBI 
//{
//   static MBCore instance;
//   return &instance;
//}
