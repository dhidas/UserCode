#ifndef GUARD_DHidasJSON_h
#define GUARD_DHidasJSON_h


#include <iostream>
#include <fstream>
#include <string>
#include <map>


class DHidasJSON
{
  public:
    DHidasJSON ();
    DHidasJSON (std::string const&, bool);
    DHidasJSON (std::string const&);
    ~DHidasJSON ();

    bool ReadFile (std::string const&, bool const Use = true);
    bool IsGoodLumiSection(int const, int const);
    void UseJSON (bool const);

  private:
    std::multimap<int, std::pair<int, int> > fMap;
    bool fUseJSON;

};


#endif
