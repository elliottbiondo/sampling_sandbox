/*
 * Contents:  Random-number sampling using the Walker-Vose alias method,
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
 * M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
 * A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)
 */
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <iostream>
//#include <sstream>
//#include <iomanip> // for setprecision
//#include <limits> // for min/max values
//#include <math.h>
//#include <vector>

#include "moab/Core.hpp"
//#include "MBTagConventions.hpp"
#include "moab/Range.hpp"


#include "sampling.h"
//MBInterface *MBI();


/*
std::vector<double> pdfFromIMesh(imesh, char* tagName)
{

}

std::vector<double> pdfFromIMesh(imesh, char* tagName, char* biasTagName)
{

}

std::vector<double> pdfFromHDF5(hdf5, char* tagName, char* biasTagName)
{

}

std::vector<double> pdfFromHDF5(hdf5, char* tagName)
{

}
*/

AliasTable::AliasTable(std::vector<double> p){
    
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

int AliasTable::drawSample(double ran1, double ran2){
    int i = (int) n * ran1;
    return ran2 < prob[i] ? i : alias[i];
}

int main()
{

  double my_array[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> my_vec(&my_array[0], &my_array[0]+5);

  AliasTable myTable(my_vec);

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

  //moab::Core *mb;
  moab::Interface *mb = new moab::Core;
  moab::ErrorCode rval;
  rval = mb->load_mesh("test.h5m");
  moab::Range hex;
  rval = mb->get_entities_by_type( 0, moab::MBHEX, hex );

  moab::Tag idxTag;
  moab::Tag phtnSrcTag;
  moab::Tag phtnSrcTag2;

  //rval = mb->tag_get_handle( "ve_idx", 0, tag_type, idxTag, moab::MB_TAG_ANY);
  rval = mb->tag_get_handle( "idx", 0, moab::MB_TYPE_INTEGER, idxTag);
  rval = mb->tag_get_handle( "phtn_src", 0, moab::MB_TYPE_DOUBLE, phtnSrcTag);
  rval = mb->tag_get_handle( "phtn_src2", 0, moab::MB_TYPE_DOUBLE, phtnSrcTag2);

  int idxTagSize;
  int phtnSrcTagSize;
  int phtnSrcTagSize2;

  rval = mb->tag_get_bytes(idxTag, idxTagSize);
  rval = mb->tag_get_bytes(phtnSrcTag, *(&phtnSrcTagSize));
  rval = mb->tag_get_bytes(phtnSrcTag2, *(&phtnSrcTagSize2));

  std::cout<< phtnSrcTagSize << std::endl;
  std::cout<< phtnSrcTagSize2 << std::endl;

  //moab::DataType tag_type;
  //rval = mb->tag_get_data_type(idxTag, tag_type);

  std::vector<int>* idxData = new std::vector<int> [idxTagSize * hex.size() / sizeof(int)];
  std::vector<double>* phtnSrcData = new std::vector<double> [phtnSrcTagSize * hex.size() / sizeof(double)];
//  std::vector<double>* phtnSrcData2 = new std::vector<double> [phtnSrcTagSize2 * (hex.size() + 2) / sizeof(double)];

  //std::vector<double>* phtnSrcData2 = new std::vector<double> [phtnSrcTagSize2 * hex.size() / 4];

  std::vector<double> phtnSrcData2;
  phtnSrcData2.resize(hex.size()*phtnSrcTagSize2/sizeof(double)); 
 //void* ptr = &phtnSrcData2[0];

//printf("\n hex.size(): %i", (int) hex.size());
//printf("\n vector length: %i double %i", (int) phtnSrcTagSize2, (int) sizeof(double));

  rval = mb->tag_get_data( idxTag, hex, idxData);
  rval = mb->tag_get_data( phtnSrcTag, hex, phtnSrcData);
  rval = mb->tag_get_data( phtnSrcTag2, hex, &phtnSrcData2[0]);

  //std::cout << phtnSrcData2.size() <<std::endl;
  
//  for(i=0; i<hex.size(); ++i){
//    printf("\n idx: %i phtn_src: %f", ((int*)(idxData))[i], ((double*)(phtnSrcData))[i]);
//  }

  for(i=0; i<hex.size()*2; ++i){
    //printf("\n phtn_src2: %f",phtnSrcData2[i]);
    std::cout << phtnSrcData2[i] << std::endl;
  }
  
  //std::vector<moab::EntityHandle> ents; 
  return 0;
}

//moab::MBInterface *MBI(){
//    static MBCore instance; 
//    return &instance;
//}

