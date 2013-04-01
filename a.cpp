#include "array.h"
#include "input-class.h"
#include <algorithm>
#include <cmath>

class C { public: int x; }; 

class CNest {
public:
  array<int> c;
  CNest() : c(6,7) {}
};

void test_array() {
  array<C> bar;
  bar.resize(5);
  bar[4].x = 12345;
  assert(bar[4].x == 12345);

  array<CNest> barNest(5,4);
  barNest(4,3).c(5,6) = -7;
  assert(barNest(4,3).c(5,6) == -7);

  array<float> foo(3,4,5);
  foo(2,3,4) = 42.0;
  // This throws "index out of bounds":  foo(4,3,2) = 44;
  assert( foo(2,3,4) == 42.0);

  float* p = (float*)foo; p[0] = foo[0];
  // This should fail to compile: int* p = (int*)foo;
}

// Build up unit tests?  In eigen, too?
// http://eigen.tuxfamily.org



// Return a's maximum amplitude.
short feature(array<short>& a) {
  const short* p = (const short*)a;
  return a.empty() ? 0 : *std::max_element(p, p+a.size());
  // short ampl = 0; for (int i=0; i<a.size(); ++i) ampl = std::max(ampl, short(std::abs(a(i)))); return ampl;
}

void plot(short x) {
  std::cout << "The feature is " << x << ".\n";
}

int main() {
  test_array();

  // Test implicit assignment op and copy constructor.
  input_t testCopying;
  input_t testAssigning(testCopying);
  input_t in = testAssigning;

  while (in) {
    array<short>& x = in.readSamples(-256, 0, 0b00000111);
    short y = feature(x);
    plot(y);
    x = in.readSeconds(0.1, 0, 0b00000111);
    y = feature(x);
    plot(y);
  }
  return 0;
}
