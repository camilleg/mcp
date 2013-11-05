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
  // operator[]
  array<C> bar;
  bar.resize(5);
  bar[4].x = 12345;
  assert(bar[4].x == 12345);

  // operator()
  array<CNest> barNest(5,4);
  barNest(4,3).c(5,6) = -7;
  assert(barNest(4,3).c(5,6) == -7);

  array<float> foo(3,4,5);
  foo(2,3,4) = 42.0;
  // This should throw "index out of bounds":  foo(4,3,2) = 44;
  assert( foo(2,3,4) == 42.0);

  float* p = (float*)foo; p[0] = foo[0];
  // This should fail to compile: int* p = (int*)foo;

  // This shold throw "negative dimension":  array<short> asdf(-256);
}

// Build up unit tests?  In eigen, too?
// http://eigen.tuxfamily.org



// Return a's maximum positive-or-negative amplitude.
short feature(array<short>& a) {
  const short* p = (const short*)a;
  if (a.empty())
    return 0;
  const std::pair<const short*,const short*> negpos = std::minmax_element(p, p+a.size());
  const short& u = *negpos.first;
  const short& v = *negpos.second;
  // Handle all cases: 0<u<v, u<0<v, u<v<0.
  return u<=0 ? v : v<=0 ? -u : std::max(-u,int(v));
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

  array<short> x(128); // intentionally too small
  while (in) {
    in.readSamples(x, -42, 0, 0x7);
    short y = feature(x);
    plot(y);
    in.readSamples(x, 256, 0, 0x3);
    y = feature(x);
    plot(y);
    in.readSeconds(x, 0.1, 0, 0x7);
    y = feature(x);
    plot(y);
  }
  return 0;
}
