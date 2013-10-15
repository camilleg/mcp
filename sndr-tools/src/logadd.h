#ifndef __LOGADD_H__
#define __LOGADD_H__

#include <algorithm> // std::max
#include <cmath>

// Log-add y to x.  For log-sum-exp idiom, to avoid underflow and overflow.
// Class T should be compatible with float or double (+, +=, <, exp, isnan, fabs).
template <class T> inline void logadd( T& x, const T& y)
{
  if (y == -INFINITY) {
    // x remains unchanged, trivially.
    return;
  }
  if (x == -INFINITY && isnan(y)) {
    // y might be -HUGE_VAL, from log(0.0).
    return;
  }

  if (isnan(x) || isnan(y))
    return;

  const T z = fabs(x-y);
  if (z > 30.0) {
    x = std::max(x, y);
    return;
  }
  if (x == y) {		// Is this special case worth testing for?
    x += log(2.0);
    return;
  }
  x = std::max(x, y) + log1p(exp(-z));
}

#endif
