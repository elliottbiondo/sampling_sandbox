/*
 * Library:   ransampl (random number sampling)
 *
 * File:      ransampl.c
 *
 * Contents:  Random-number sampling using the Walker-Vose alias method,
 *            as described by Keith Schwarz (2011)
 *            [http://www.keithschwarz.com/darts-dice-coins]
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
 *
 * License:   see ../COPYING (FreeBSD)
 * 
 * Homepage:  apps.jcns.fz-juelich.de/ransampl
 */

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <iostream>
#include <sstream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <math.h>
#include <vector>

#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Core.hpp"

//MBInterface *MBI();


/*
std::vector<double> pdfFromIMesh(imesh, char* tagName){

}

std::vector<double> pdfFromIMesh(imesh, char* tagName, char* biasTagName){

}

std::vector<double> pdfFromHDF5(hdf5, char* tagName, char* biasTagName){

}

std::vector<double> pdfFromHDF5(hdf5, char* tagName){

}
*/

class AliasTable{
  std::vector<double> prob;
  std::vector<int> alias;
  int n;
public:
  AliasTable(std::vector<double> p);
  int drawSample(double ran1, double ran2);
};

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

int main(){

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
    printf("%i    %f   %f \n", i+1, (float) answers[i]/N, (float) (i+1)/15.0);
  }

  moab::Core *mb;
  moab::ErrorCode rval = mb->load_mesh("test.h5m");
  moab::Tag phtnSrcTag;

  return 0;
}

//MBInterface *MBI(){
//    static MBCore instance;
//    return &instance;
//}
