#include "DHidasLJAna/LeptonPlusJets/interface/DHidasJSON.h"


DHidasJSON::DHidasJSON ()
{
}


DHidasJSON::DHidasJSON (std::string const& InFileName)
{
  ReadFile(InFileName);
}


DHidasJSON::~DHidasJSON ()
{
}


bool DHidasJSON::ReadFile (std::string const& InFileName)
{
  std::ifstream f(InFileName.c_str());
  if (!f.is_open()) {
    throw;
  }

  char c;
  int n;

  int run, lb, le;
  bool startlb = false;
  while(!f.eof()) {
    c = f.peek();
    while(!f.eof() && !( (c >= '0') && (c <= '9'))) {
      c = f.get();
      c = f.peek();
    }
    f >> n;
    if(n > 100000) {
      run = n;
    } else {
      if(!startlb) {
        lb = n;
        startlb = true;
      } else {
        le = n;
        startlb = false;
        fMap.insert( std::make_pair<int, std::pair<int, int> >(run, std::make_pair<int, int>(lb, le) ) );
      }
    }
  }
  std::cout << "Number of good lumi sections: " << fMap.size() << std::endl;

  return true;
}


bool DHidasJSON::IsGoodLumiSection(int const run, int const lumis)
{
  std::multimap<int, std::pair<int, int> >::iterator p = fMap.find(run);
  if (p == fMap.end()) {
    return false;
  }

  for (std::multimap<int, std::pair<int, int> >::iterator b = fMap.upper_bound(run); p != b; ++p) {
    if (lumis >= p->second.first && lumis <= p->second.second) {
      return true;
    }
  }
  return false;
}

