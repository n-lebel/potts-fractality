#include "key.h"

Key::Key(float x, float y): i(x), j(y) {}

bool Key::operator<(const Key &other) const {
      if (i < other.i) return true;
      if (other.i < i) return false;
      return j < other.j;
}

bool Key::operator==(const Key &other) const {
	if (i == other.i && j == other.j) return true;
	else return false;
}

bool Key::operator!=(const Key &other) const {
	if (i == other.i && j == other.j) return false;
	else return true;
}

std::ostream& operator<<(std::ostream& os, const Key& k) {

	os << "(" << k.i << ", " << k.j << ")";
	return os;
}
