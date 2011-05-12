#ifndef GUARD_DHRunTracker_h
#define GUARD_DHRunTracker_h

#include <vector>
#include <algorithm>

class DHRunTracker
{
  public:
    DHRunTracker ();
    ~DHRunTracker ();

    bool IsDuplicate(unsigned int);

  private:
    std::vector<unsigned int> fVector;
};


#endif
