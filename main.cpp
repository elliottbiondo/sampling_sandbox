#include "sampling.hpp"
#include <iostream>

int main(int argc, char* argv[]){

  int i, j;
  /* Working alias table if you make it public again.
  double my_array[5] = {0.066667, 0.133333, 0.200000, 0.266667, 0.333333};
  std::vector<double> my_vec(&my_array[0], &my_array[0]+5);
  Sampling::AliasTable myTable(my_vec);
  int answers[] = {0, 0, 0, 0, 0};
  int N = 1000000;
  double rand1, rand2;
  for(i=0; i<N; i++){
     rand1 = (double) rand()/RAND_MAX;
     rand2 = (double) rand()/RAND_MAX;
     answers[myTable.sample_pdf(rand1, rand2]++;
  }
  printf("bin |  prob  | expected prob\n");
  for(i=0; i<5; i++){
    printf("%i    %f   %f \n", i+1, (double) answers[i]/N, (double) (i+1)/15.0);
  }
 */


 //Sampling& sampling = *Sampling::instance();
 gggsampling_setup_();
 //fsampling_setup_(argv[1], argv[2], argv[3], false);
 //sampling.sampling_setup_(argv[1], argv[2], argv[3], true, argv[4]);

 double rands[6];
 double x, y, z, e, w;
 //sampling.particle_birth(rands, x, y, z, E);
  std::ofstream myfile;
  myfile.open ("samples.out");
 for(i=0; i<5000; i++){
   for(j=0; j<6; j++){
     rands[j] = (double) rand()/RAND_MAX;
   }
   fparticle_birth_(rands, x, y, z, e, w);
   //myfile << x <<" "<< y <<" "<<z<< " "<< e << " "<< w <<std::endl;
 }
return 0;
}
