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

Sampling::AliasTable::AliasTable(std::vector<double> p){
    n = p.size();
    prob.resize(n);
    alias.resize(n);
    std::vector<double> S(n);
    std::vector<double> L(n);
    int i, a, g;

    for(i=0; i<n; ++i) 
      p[i] *= n;

    // Set separate index lists for small and large probabilities:
    int nS = 0;
    int nL = 0;
    for(i=n-1; i>=0; --i) {
        // at variance from Schwarz, we revert the index order
        if(p[i] < 1)
            S[nS++] = i;
        else
            L[nL++] = i;
    }

    // Work through index lists
    while(nS && nL){
        a = S[--nS]; // Schwarz's l
        g = L[--nL]; // Schwarz's g
        prob[a] = p[a];
        alias[a] = g;
        p[g] = p[g] + p[a] - 1;
        if (p[g] < 1)
            S[nS++] = g;
        else
            L[nL++] = g;
    }

    while(nL)
        prob[L[--nL]] = 1;

    while(nS)
        // can only happen through numeric instability
        prob[ S[--nS] ] = 1;
}

int Sampling::AliasTable::draw_sample(double rand1, double rand2){
    int i = (int) n * rand1;
    return rand2 < prob[i] ? i : alias[i];
}

void Sampling::SamplingSetup(char* fileName, char* src_tag_name, char* bias_tag_name, char* e_bounds_tag_name){


}




void Sampling::SamplingSetup(char* fileName, char* src_tag_name, char* e_bounds_tag_name){
  MBErrorCode rval;
  MBEntityHandle loaded_file_set;
  rval = MBI->create_meshset(MESHSET_SET, loaded_file_set );
  //assert( rval == MB_SUCCESS );
  rval = MBI->load_file( fileName, &loaded_file_set );
  //assert( rval == MB_SUCCESS );
  rval = MBI->get_entities_by_dimension(loaded_file_set, 3, ves);

  int num_hex, num_tet;
  rval = MBI->get_number_entities_by_type(loaded_file_set, MBHEX, num_hex);
  rval = MBI->get_number_entities_by_type(loaded_file_set, MBTET, num_tet);
  if(num_hex == ves.size()){
    ve_type = MBHEX;
    verts_per_vol = 8;
  } else if (num_tet == ves.size()){
    ve_type = MBTET;
    verts_per_vol = 4;
  }
  else exit(1);

  std::vector<double> volumes =  find_volumes();

  //assert( rval == MB_SUCCESS );
  // get tag handle
  MBTag src_tag, bias_tag;
  rval = MBI->tag_get_handle(src_tag_name, moab::MB_TAG_VARLEN, MB_TYPE_DOUBLE, src_tag);
  // THIS ASSERT FAILS because we do not know number of energy groups a priori.
  //assert( rval == MB_SUCCESS );
  // get tag size
  int tag_size;
  rval = MBI->tag_get_bytes(src_tag, *(&tag_size));
  //assert( rval == MB_SUCCESS );
  tag_len = tag_size/sizeof(double);

  std::vector<double> pdf(ves.size()*tag_len); 
  rval = MBI->tag_get_data(src_tag, ves, &pdf[0]);
  //assert( rval == MB_SUCCESS );

  /*
  MBTag e_tag;
  std::cout << e_bounds_tag_name << std::endl;
  rval = MBI->tag_get_handle(e_bounds_tag_name, 3, MB_TYPE_DOUBLE, e_tag);
  std::cout << rval << std::endl;
  assert(rval == MB_SUCCESS);
  rval = MBI->tag_get_data(e_tag, &loaded_file_set, 1, &e_bounds[0]);
  std::cout << rval << std::endl;
  */
  //assert(rval == MB_SUCCESS);
//  std::cout << e_bounds[0] << std::endl;
  //std::cout << e_bounds[0] << e_bounds[1] << e_bounds[2] << std::endl;

  int i, j;
  for(i=0; i<ves.size(); ++i){
    for(j=0; j<tag_len; ++j){
     pdf[i*tag_len + j] *=  volumes[i];
    }
  }
  at = new AliasTable(pdf);
}

std::vector<double> Sampling::find_volumes(){
  std::vector<double> volumes (ves.size());
  MBErrorCode rval;
  std::vector<MBEntityHandle> connect;
  rval = MBI->get_connectivity_by_type(ve_type, connect);
  double coords[verts_per_vol*3];
  int i;
  for(i=0; i<ves.size(); ++i){
    rval=MBI->get_coords(&connect[verts_per_vol*i], verts_per_vol, &coords[0]);
    volumes[i] = measure(ve_type, verts_per_vol, &coords[0]);

    if(ve_type == MBHEX){
      MBCartVect o(coords[0], coords[1], coords[2]);
      MBCartVect x(coords[3], coords[4], coords[5]);
      MBCartVect y(coords[9], coords[10], coords[11]);
      MBCartVect z(coords[12], coords[13], coords[14]);
      vector_points vp = {o, x-o, y-o, z-o};
      cart_sampler.push_back(vp);
   }else if (ve_type == MBTET){
      MBCartVect o(coords[0], coords[1], coords[2]);
      MBCartVect x(coords[3], coords[4], coords[5]);
      MBCartVect y(coords[6], coords[7], coords[8]);
      MBCartVect z(coords[9], coords[10], coords[11]);
      vector_points vp = {o, x-o, y-o, z-o};
      cart_sampler.push_back(vp);
    }
  }
  return volumes;
}


void Sampling::SampleXYZE(double* rands, double &x, double &y, double &z, double &e, double &w){
  // get indices
  int pdf_idx = at->draw_sample(rands[0], rands[1]);
  int ve_idx = pdf_idx/tag_len;
  int e_idx = pdf_idx % tag_len;
  
  // get x, y, z
  if(ve_type == MBHEX){
    get_xyz(ve_idx, &rands[2],x,y,z);
  } else if (ve_type == MBTET){
    double s = rands[2];
    double t = rands[3];
    double u = rands[4];

    if(s + t > 1){
      s = 1.0 - s;
      t = 1.0 - t;
    }
    
    if(s + t + u > 1){
      if(t + u > 1){
        double old_t = u;
        t = 1.0 - u;
        u = 1.0 - s - old_t;
      }else if (t + u <= 1){
        double old_s = s;
        s = 1.0 - t - u;
        u = old_s + t + u - 1;
      }
    }
    double new_rands[3] = {s, t, u};
    get_xyz(ve_idx, new_rands, x, y, z);
  }

  get_e(e_idx, &rands[5], e);
  w = 1.0;
}

void Sampling::get_xyz(int ve_idx, double* rands, double &x, double &y, double &z){

  MBCartVect a = rands[0]*cart_sampler[ve_idx].x_vec + \
                 rands[1]*cart_sampler[ve_idx].y_vec + \
                 rands[2]*cart_sampler[ve_idx].z_vec + \
                 cart_sampler[ve_idx].o_point;
  x = a[0];
  y = a[1];
  z = a[2];
}

void Sampling::get_e(int e_idx, double* rand, double &e){
   e_bounds.push_back(1.1);
   e_bounds.push_back(1.2);
   e_bounds.push_back(1.3);
   double e_min = e_bounds[e_idx];
   double e_max = e_bounds[e_idx + 1];
   e = rand[0] * (e_max - e_min) + e_min;

}

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
     answers[myTable.draw_sample(rand1, rand2)]++;
  }
  printf("bin |  prob  | expected prob\n");
  for(i=0; i<5; i++){
    printf("%i    %f   %f \n", i+1, (double) answers[i]/N, (double) (i+1)/15.0);
  }
 */




 Sampling& sampling = *Sampling::instance();
 sampling.SamplingSetup(argv[1], argv[2], argv[3]);

 double rands[6];
 double x, y, z, e, w;
 //sampling.SampleXYZE(rands, x, y, z, E);
  std::ofstream myfile;
  myfile.open ("unstr.out");
 for(i=0; i<5000; i++){
   for(j=0; j<6; j++){
     rands[j] = (double) rand()/RAND_MAX;
   }
   sampling.SampleXYZE(rands, x, y, z, e, w);
   myfile << x <<" "<< y <<" "<<z<< " "<< e << " "<< w <<std::endl;
 }

return 0;
}
