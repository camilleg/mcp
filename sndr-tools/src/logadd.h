#ifndef __LOGADD_H__
#define __LOGADD_H__

#include <algorithm> // std::max
#include <cmath>

#include "real_t.h"

// Log-add y to x.  For log-sum-exp idiom, to avoid underflow and overflow.

inline void logadd(real_t& x, const real_t& y)
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

  const real_t z = fabs(x-y);
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
