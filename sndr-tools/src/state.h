// Base class for ANN, GMM, etc.
// language: C++
// author  : Camille Goudeseune

#ifndef __STATE_H__
#define __STATE_H__

#include "real_t.h"

class state_t {
public:
  // Learn data "x".
  virtual void train(const array<real_t> &x, const int iters = 100, const state_t &S = state_t(), bool prior = false) {}

  virtual void save(const std::string& filename) const {}
  virtual void load(const std::string& filename) {}

private:
  // Evaluate log likelihoods of M*N data x into N*K-array p.
  virtual void likelihoods(array<real_t> &p, const array<real_t> &x) {}
};

#endif
