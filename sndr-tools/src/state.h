// "Base" class for ANN, GMM, etc.
// language: C++
// author  : Camille Goudeseune

#ifndef __STATE_H__
#define __STATE_H__

template <class T>
class state_t {
  // These methods aren't virtual, because a template class may not have virtual functions.
  // Instead, classifier subclasses (gmm_t, ann_t) just override these (empty) methods with their own.
public:
  // Learn data "x".
  void train(const array<T> &x, const int iters = 100, const state_t<T> &S = state_t<T>(), bool prior = false) {}

  void save(const std::string& filename) const {}
  void load(const std::string& filename) {}

private:
  // Evaluate log likelihoods of M*N data x into N*K-array p.
  void likelihoods(array<T> &p, const array<T> &x) {}
};

#endif
