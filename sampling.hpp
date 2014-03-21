#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "DagMC.hpp"
#include "MBTagConventions.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"
#include "moab/Core.hpp"
#include "moab/GeomUtil.hpp"
#include "MBCore.hpp"
#include "measure.hpp"
#include "MBCartVect.hpp"

class Sampling
{
public:
  static Sampling *instance(MBInterface *mb_impl = NULL);
  ~Sampling();
  void SamplingSetup(char* fileName, char* tagName);
  void SampleXYZE(const double* rands, double &x, double &y, double &z, double &E);

  class AliasTable
  {
  private:
    std::vector<double> prob;
    std::vector<int> alias;
    int n;

    friend class Sampling;
  
  public:
    int draw_sample(double ran1, double ran2);
    AliasTable(std::vector<double> p);
  };

private:
  //functions
  void pdfFromMesh(char* fileName, char* tagName);
  std::vector<double> find_volumes();
  //variable
  int tag_len;
  MBRange ves;
  MBErrorCode rval;
  std::vector<double> pdf;
  MBEntityType ve_type;
  int verts_per_vol;
  AliasTable* at;
  struct vector_points{
    MBCartVect o_point;
    MBCartVect x_vec;
    MBCartVect y_vec;
    MBCartVect z_vec;
  };
  std::vector<vector_points> cart_sampler;

public:
  MBInterface* moab_instance() {return mbImpl;}

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
