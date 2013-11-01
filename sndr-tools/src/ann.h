// ann.h -- Artificial Neural Network class
// language: C++
// author  : Camille Goudeseune

#ifndef __ANN_H__
#define __ANN_H__

#include "array.h"
#include "state.h"

#include <fstream>
#include <iostream>

template <class T>
class annLayer_t {
  void g(array<T>& dst, const array<T>& src) {}      // Layer's nonlinearity.  dst and src are vectors.
  void gprime(array<T>& dst, const array<T>& src) {} // Derivative of g().
};
// Subclasses: tanhLayer logisticLayer, softmaxLayer, rbfLayer, linearLayer, perceptronLayer, hingeLayer, stepLayer
// These override g() and gprime().

template <class T>
class annMetric_t {
  void e(array<T>& dst, const array<T>& t, const array<T>& y) {} // Error of network outputs y with respect to targets t.
  void eprime(array<T>& dst, const array<T>& t, const array<T>& y) {} // Derivative of e().  Delta k.

};
// Subclasses: l2Metric, lpMetric, softmaxEntropyMetric, binaryEntropyMetric.
// These override e() and eprime().


// Gaussian mixture model class
class ann_t: public state_t {
public:

  // Learn data "x"
  void train( const array<real_t> &x, const int iters = 100, const ann_t &G = ann_t(), bool prior = false)
  {
  }

  void save(const std::string& filename) const
  {
  }

  void load(const std::string& filename)
  {
  }

  // Compute output of every layer (outputs is 2-D).
  void forwardPropagate(array<real_t>& outputs, const array<real_t>& inputToFirstLayer)
  {
    // For each member of layers, call its ->g().
  }

  // Compute gradient of error of tokenVector w.r.t all weightMatrices.
  void gradient(array<real_t>& outputs, const array<real_t>& tokenVector, const array<real_t>& x, const array<real_t>& Z)
  {
    // During back propagation, call metric.eprime() and layers[].gprime().
  }

private:
  // Evaluate log likelihoods of M*N data x into N*K-array p.
  void likelihoods(array<real_t> &p, const array<real_t> &x)
  {
  }

  array<annLayer_t<real_t> > layers;
  array<real_t> weightMatrices; // one per layer of the NN
  annMetric_t<real_t> metric;
};

#endif
