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

class Sampling
{
public:
  static Sampling *instance(MBInterface *mb_impl = NULL);
  ~Sampling();
  //functions
  //std::vector<double> pdfFromMesh(char* fileName, char* tagName);
  void pdfFromMesh(char* fileName, char* tagName);
  std::vector<double> find_volumes();
  //variable
  
  MBTag phtnSrcTag;
  int tagLen;
  MBTag idxTag;
  MBRange ves;
  MBErrorCode rval;
  std::vector<double> phtnSrcData;
  MBEntityType ve_type;
  int verts_per_vol;

public:
  class AliasTable
  {
  private:
    std::vector<double> prob;
    std::vector<int> alias;
    int n;

    friend class Sampling;
  
  public:
    int drawSample(double ran1, double ran2);
    AliasTable(std::vector<double> p);
  };

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
