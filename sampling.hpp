
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>


#include "MBCore.hpp"

class Sampling
{
public:
  class AliasTable
  {
    std::vector<double> prob;
    std::vector<int> alias;
    int n;
  
  public:
    AliasTable(std::vector<double> p);
    int drawSample(double ran1, double ran2);
  };
  
};
