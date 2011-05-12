#ifndef GUARD_DHRunTracker_h
#define GUARD_DHRunTracker_h

#include <vector>
#include <set>
#include <map>
#include <algorithm>

class DHRunTracker
{
  public:
    DHRunTracker ();
    ~DHRunTracker ();

    bool IsDuplicate(unsigned int, unsigned int);

  private:
    std::set< std::pair<unsigned int, unsigned int> > fSet;
};


#endif
