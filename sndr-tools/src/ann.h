// ann.h -- Artificial Neural Network class
// language: C++
// author  : Camille Goudeseune

#ifndef __ANN_H__
#define __ANN_H__

#include <fstream>
#include <iostream>

#include "array.h"
#include "state.h"

// Gaussian mixture model class
template <class T>
class ann_t: public state_t<T> {
public:

  // Learn data "x"
  void train( const array<T> &x, const int iters = 100, const ann_t<T> &G = ann_t<T>(), bool prior = false)
  {
  }

  void save(const std::string& filename) const
  {
  }

  void load(const std::string& filename)
  {
  }

private:
  // Evaluate log likelihoods of M*N data x into N*K-array p.
  void likelihoods(array<T> &p, const array<T> &x)
  {
  }
};

#endif
