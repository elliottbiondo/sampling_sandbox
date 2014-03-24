#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "moab/Range.hpp"
#include "MBCore.hpp"
#include "measure.hpp"
#include "MBCartVect.hpp"

class Sampling
{
public:
  static Sampling *instance(MBInterface *mb_impl = NULL);
  MBInterface* moab_instance() {return mbImpl;}
  ~Sampling();
  void SamplingSetup(char* fileName, char* src_tag_name, char* e_bound_tag_name);
  void SamplingSetup(char* fileName, char* src_tag_name, char* e_bounds_tag_name, char* bias_tag_name);
  void SampleXYZE(double* rands, double &x, double &y, double &z, double &e, double &w);


private:
  //functions
  std::vector<double> find_volumes();
  //variable
  bool bias;
  int tag_len;
  MBRange ves;
  MBEntityType ve_type;
  std::vector<double> e_bounds;
  int verts_per_vol;
  struct vector_points{
    MBCartVect o_point;
    MBCartVect x_vec;
    MBCartVect y_vec;
    MBCartVect z_vec;
  };
  std::vector<vector_points> cart_sampler;
  void get_xyz(int ve_idx, double* rands, double &x, double &y, double &z);
  void get_e(int e_idx, double* rand, double &e);
  class AliasTable
  {
  private:
    std::vector<double> prob;
    std::vector<int> alias;
    int n;
  public:
    int draw_sample(double ran1, double ran2);
    AliasTable(std::vector<double> p);
  };
  AliasTable* at;

private:
  Sampling(MBInterface *mb_impl);
  static void create_instance(MBInterface *mb_impl = NULL);
  static Sampling *instance_;
  MBInterface *mbImpl;
};


inline Sampling *Sampling::instance(MBInterface *mb_impl)
{
  if (NULL == instance_) create_instance(mb_impl);
  return instance_;
}
