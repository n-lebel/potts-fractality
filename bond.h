#ifndef bond_h
#define bond_h

#include "key.h"

class Bond {
public:
   Key a;
   Key b;

   Bond(Key one, Key two);
   bool operator<(const Bond &other) const;
   bool operator==(const Bond &other) const;
};
std::ostream& operator<<(std::ostream& os, const Bond& b);

#endif
