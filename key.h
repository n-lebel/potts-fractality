#ifndef key_h
#define key_h

#include <map>
#include <fstream>

class Key {
public:
   float i;
   float j;

   Key(float x, float y);
   bool operator<(const Key &other) const;
   bool operator==(const Key &other) const;
   bool operator!=(const Key &other) const;
};

std::ostream& operator<<(std::ostream& os, const Key& k);
#endif
