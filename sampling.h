class AliasTable
{
  std::vector<double> prob;
  std::vector<int> alias;
  int n;

public:
  AliasTable(std::vector<double> p);
  int drawSample(double ran1, double ran2);
};

