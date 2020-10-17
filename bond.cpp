#include "bond.h"

Bond::Bond(Key one, Key two): a(one), b(two) {}

bool Bond::operator<(const Bond &other) const {
      if (a < other.a) return true;
      if (other.a < a) return false;
      return b < other.b;


}

bool Bond::operator==(const Bond &other) const {
	if (a == other.a && b == other.b) return true;
	if (b == other.a && a == other.b) return true;
	else return false;

}

std::ostream& operator<<(std::ostream& os, const Bond& l) {
	os << "(" << l.a << "--" << l.b << ")";
	return os;
}
