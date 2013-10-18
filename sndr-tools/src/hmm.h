// hmm.h -- Hidden Markov Model class with Gaussian Mixture Model states
// language: C++
// author  : Paris Smaragdis

#ifndef __HMM_H__
#define __HMM_H__

#include <fstream>
#include <iostream>

#include "ann.h"
#include "gmm.h"
#include "array.h"
#include "logadd.h"

// Hidden Markov Model class

template <class T>
class hmm_t {
public:
  int S; // States
  int K; // Gaussians per state

  array<T> lPi; // Model parameter: initial probabilities
  array<T> lA;  // Model parameter: matrix of transition probabilities

  // State model parameters.  Could be a GMM, HMM, ANN, etc.
  array<gmm_t<T> > smps;
  array<ann_t<T> > smps_unused_but_at_least_it_compiles;

private:
  // For each element of smps.
  class statemodel_t {
  public:
    array<T> g;      // (From line 151 in sndr-tools_v0.6/src/hmm.h.  A "buffer.")
    array<T> la, lb; // Alpha and beta parameters for iterations in Baum-Welch training.
  };
  array<statemodel_t> statemodel; // Related stuff for each element of smps
  array<T> xi, txi; // SxS parameters for Baum-Welch iterations

public:
  // Constructor
  hmm_t(const int s=0, const int k=1): S(s), K(k) {}

  // Learn data
  void train(const array<T> &x, const int iters = 100, const hmm_t<T> &H = hmm_t<T>())
  {
    // Reverse x's dims.
    const int M = x.n;
    const int N = x.m;
    std::cout << "HMM training " << M << " x " << N << std::endl;

    // Setup state model parameters
    smps.resize(S);
    for (int s=0; s<S; ++s) {
      // M and N come from arg x.  K comes from constructor, or from load().
      smps(s).init(M,K);
    }

    if (H.S == 0) {
      // Initial values of Gaussians
      for (int s=0; s<S; ++s)
	smps(s).init2();

      // Initial values of initial and transition probabilities
      lPi.resize(S);
      for (int s=0; s<S; ++s) lPi(s) = log(1.0/S);
      lA.resize(S, S);
      for (int i=0; i<S*S; ++i) lA[i] = log(T(rand())/RAND_MAX);
      for (int i=0; i<S; ++i) {
	T ls = log(0.0);
	for (int j=0; j<S; ++j) logadd(ls, lA(i,j));
	for (int j=0; j<S; ++j) lA(i,j) -= ls;
      }
    }else{
      // Copy values of Gaussians
      smps = H.smps;
      for (int s=0; s<S; ++s)
	if (smps(s).uninitialized())
	  throw std::runtime_error( "hmm_t::train(): uninitialized gmm.");

      // Copy values of initial and transition probabilities
      lPi = H.lPi;
      lA = H.lA;
    }

    // Allocate buffers
    array<T> q( N, K, S);
    array<T> lp( S, N);
    array<T> lk( iters);
    statemodel.resize(S);
    for( int s=0; s<S; ++s) {
      statemodel_t& store = statemodel(s); // related stuff for each element of smps
      store.g.resize(N,K);
      store.la.resize(N);
      store.lb.resize(N);
    }
    xi .resize(S, S);
    txi.resize(S, S);

    array<int> dummy(K);
    for (int i=0; i<K; ++i)
      dummy[i] = 1;

    // Iterate expectation-maximization
    for( int it=0; it<iters; ++it){

      // *** Expectation step ***

      // Get likelihoods from each gaussian from each state
#pragma omp parallel for
      for( int s=0; s<S; ++s)
	for( int k=0; k<K; ++k)
	  for( int j=0; j<N; ++j)
	    q(j,k,s) = smps(s).placeholderForFunctionName(x,j,k);

      // Compute overall state likelihoods
      for( int s=0; s<S; ++s)
	for( int j=0; j<N; ++j){
	  T tp = log( 0.0);
	  for( int k=0; k<K; ++k)
	  logadd( tp, q(j,k,s));
	  lp(s,j) = tp;
	}

      // Get alphas
      for (int s=0; s<S; ++s)
	statemodel(s).la(0) = lPi(s) + lp(s,0);
      for (int t=0; t<N-1; ++t) {
	for (int j=0; j<S; ++j) {
	  T ls = log(0.0);
	  for (int s=0; s<S; ++s)
	    logadd( ls, statemodel(s).la(t) + lA(s,j));
	  statemodel(j).la(t+1) = ls + lp(j,t+1);
	}
      }

      // Get betas
      for (int s=0; s<S-1; ++s)
	statemodel(s).lb(N-1) = log(0.0);
      statemodel(S-1).lb(N-1) = 0;
      for (int t = N-2; t >= 0; --t) {
	for (int i=0; i<S; ++i) {
	  T ls = log(0.0);
	  for (int j=0; j<S; ++j)
	    logadd(ls, statemodel(j).lb(t+1) + lA(i,j) + lp(j,t+1));
	  statemodel(i).lb(t) = ls;
	}
      }

      // Get Xi
      for( int i=0; i<S*S; ++i)
	xi[i] = log(0.0);
      for (int t=0; t<N-1; ++t) {
	T ls = log(0.0);
	for (int s=0; s<S; ++s) {
	  for (int j=0; j<S; ++j) {
	    txi(s,j) = lA(s,j) + statemodel(s).la(t) + lp(j,t+1) + statemodel(j).lb(t+1);
	    logadd(ls, txi(s,j));
	  }
	}
	for( int i=0; i<S*S; ++i)
	  logadd(xi[i], txi[i] - ls);
      }

      // Get gamma
      array<T> la_lb(S);
      for( int j=0; j<N; ++j) {
	T ls = log(0.0);
	for (int s=0; s<S; ++s) {
	  la_lb[s] = statemodel(s).la(j) + statemodel(s).lb(j);
	  logadd( ls, la_lb[s]);
	}
	for (int s=0; s<S; ++s){
	  const T tg = la_lb[s] - lp(s,j) - ls;
	  for( int k=0; k<K; ++k) {
	    statemodel(s).g(j,k) = tg + q(j,k,s);
	  }
	}
      }

      // Get overall likelihood
      lk(it) = log(0.0);
      for (int s=0; s<S; ++s)
	logadd(lk(it), statemodel(s).la(N-1));
      if( ((it+1)%5==0) || it == iters-1)
	std::cout << "HMM iteration " << it+1 << " of " << iters << ": likelihood " << lk(it) << std::endl;

      // Exit log domain
      for (int s=0; s<S; ++s)
	for( int i=0; i<N*K; ++i)
	  statemodel(s).g[i] = exp(statemodel(s).g[i]);


      // *** Maximization step ***

      for (int s=0; s<S; ++s) {
	// Initial probabilities
	T tp = log(0.0);
	for( int i=0; i<K; ++i)
	  logadd(tp, log(statemodel(s).g(0,i)));
	lPi(s) = tp;

	// Transition matrix
	T ls = log(0.0);
	for( int j=0; j<S; ++j)
	  logadd(ls, xi(s,j));
	for( int j=0; j<S; ++j)
	  lA(s,j) = xi(s,j) - ls;
      }

      for (int s=0; s<S; ++s)
	smps(s).maximize(statemodel(s).g, x, dummy);
    }
  }

  // Classify data "x" into "q" using a known HMM
  void classify( array<int> &q, const array<T> &x, const array<T> &bias = array<T>(), const int ist = -1)
  {
    // Input dims of array are reversed
    const int M = x.n, N = x.m;

    // Get state probabilities
    array<T> lB( S, N);
    for( int i=0; i<S*N; ++i)
      lB[i] = log(0.0); // Init all lB(s,j)'s.

#pragma omp parallel for
    for( int s=0; s<S; ++s)
      for( int k=0; k<K; ++k)
	for( int j=0; j<N; ++j)
	  logadd(lB(s,j), smps(s).placeholderForFunctionName2(x,j,k,K,M));

    // Add the bias to first state
    if (!bias.empty()) {
      if (int(bias.size()) != S)
	throw std::runtime_error( "classify(): bias vector has mismatched size");
      for (int i=0; i<S; ++i)
	for (int j=0; j<N; ++j)
	  lB(i,j) += bias(i);
    }

    viterbi( q, lB, ist);
  }

private:
  T sq(const T x) const { return x*x; }

  // Viterbi decoding
  void viterbi(array<int> &q, const array<T> &lB, const int ist = -1)
  {
    const int N = lB.n;
    array<T> d( S, 2);
    array<int> p( S, N);

    // Override initial state if one is provided
    array<double> nlPi(S);
    if( ist != -1){
      for( int i=0; i<S; ++i)
	nlPi(i) = log( 0.0);
      nlPi(ist) = 0.0;
    }else
      for( int i=0; i<S; ++i)
	nlPi(i) = lPi(i);

    // Initialize
    bool di = 0;
    for( int i=0; i<S; ++i){
      d(i,di) = nlPi(i) + lB(i,0);
      p(i,0) = 0;
    }

    // Propagate
    for( int t=1; t<N; ++t){
      for( int i=0; i<S; ++i){

	// Find max likelihood and its index on temp variables
	T dt = -HUGE_VAL;
	int pt = 0;
	for( int j=0; j<S; ++j){
	  T tt = lA(j,i) + d(j,di);
	  if( tt > dt){
	    dt = tt;
	    pt = j;
	  }
	}

	// Assign to arrays
	d(i,!di) = dt + lB(i,t);
	p(i,t) = pt;
      }
      di = !di;
    }

    // Terminate
    T l = d(0,di);
    q.resize( N);
    q(N-1) = 0;
    for( int i=1; i<S; ++i)
      if( d(i,di) > l){
	l = d(i,di);
	q(N-1) = i;
      }

    // Backtrack
    for( int t = N-2; t >= 0; --t)
      q(t) = p(q(t+1),t+1);
  }

#ifdef UNUSED
  // Short-time Viterbi decoding
  void stviterbi(array<int> &q, const array<T> &lB)
  {
    const int M = lB.m;
    const int N = lB.n;
    q.resize( N);

    // Initial probability vector
    array<T> pp(lPi);

    int a = 0;
    for( int b=1; b<N; ){

      // Get state paths for a short segment
      array<int> s;
      iviterbi( lB, a, b, pp, s);

      // Find fusion point, where all possible solutions converge to the same state
      int t = -1;
      for( int j=0; j<s.n; ++j) {
	for( int i=0; i<S-1; ++i) {
	  if (s(i,j) != s(i+1,j)) {
	    t = j;
	    j = s.n; // break out of outer loop too
	    break;
	  }
	}
      }

      if( t >= 0){
	// Keep advance
	for( int i=0; i<t; ++i)
	  q(a+i) = s(0,i);
	for( int i=0; i<S; ++i)
	  pp(i) = lA(q(a+t-1),i);
	a += t;
	b = a+1;
      }else{
	// Increase window and repeat
	++b;
      }
    }
  }

  // Interim Viterbi decoding
  void iviterbi(array<int> &q, const array<T> &lB, const array<T> &lP, const int a, const int b)
  {
    const int N = b-a;

    array<T> d( S, 2);
    array<int> p( S, N);

    // Initialize
    bool di = 0;
    for( int i=0; i<S; ++i){
      d(i,di) = lP(i) + lB(i,a);
      p(i,0) = 0;
    }

    // Propagate
    for( int t=1; t<N; ++t){
      for( int i=0; i<S; ++i){

	// Find max likelihood and its index on temp variables
	T dt = -HUGE_VAL;
	int pt = 0;
	for( int j=0; j<S; ++j){
	  T tt = lA(j,i) + d(j,di);
	  if( tt > dt){
	    dt = tt;
	    pt = j;
	  }
	}

	// Assign to arrays
	d(i,!di) = dt + lB(i,a+t);
	p(i,t) = pt;
      }
      di = !di;
    }

    // Try all possible terminations
    q.resize( S, N);
    for( int i=0; i<S; ++i){
      q(i,N-1) = i;

      // Backtrack
      for( int t = N-2; t >= 0; --t)
	q(i,t) = p(q(i,t+1),t+1);
    }
  }
#endif

public:
  void save(const std::string& filename) const
  {
    using namespace std;
    if (filename.empty())
      throw runtime_error( "hmm_t::save(\"\") failed.");
    ofstream f( filename.c_str(), ios::out | ios::binary);
    if (!f)
      throw runtime_error( "hmm_t::save('" + filename + "') failed.");

    f.write((const char*)&S, sizeof(int)); // number of states
    lPi.write(f);                          // initial log probabilities
    lA.write(f);                           // log transition matrix
    f.write((const char*)&K, sizeof(int)); // number of gaussians

    const int M = smps(0).dimensions();
    std::cout << ( "DEBUG hmm_t::save('" + filename + "'): number of dimensions M is " + to_str(M) + ".\n");
    f.write((char*)&M, sizeof(int)); // dimension
    for (int s=0; s<S; ++s)
      smps(s).writeContents(f);
    if (!f)
      throw runtime_error( "hmm_t::save('" + filename + "') failed.");
    cout << "Saved HMM file " << filename << ".\n";
  }

  void load(const std::string& filename)
  {
    using namespace std;
    if (filename.empty())
      throw runtime_error( "hmm_t::load(\"\") failed.");

    ifstream f( filename.c_str(), ios::in | ios::binary);
    if (!f)
      throw runtime_error( "hmm_t::load('" + filename + "') failed.");

    // number of states
    f.read( (char*)&S, sizeof( int));
    if (S <= 0)
      throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of states " + to_str(S) + ".");

    // initial log probabilities
    lPi.read(f, S);

    // log transition matrix
    lA.read(f, S, S);

    // number of gaussians
    f.read( (char*)&K, sizeof( int));
    if (K <= 0)
      throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of gaussians per state " + to_str(K) + ".");

    // dimension
    int M;
    f.read( (char*)&M, sizeof( int));
    if (M <= 0)
      throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of dimensions " + to_str(M) + ".");

    smps.resize(S);
    for (int s=0; s<S; ++s) smps(s).init(M, K, f);

    if (!f)
      throw runtime_error( "hmm_t::load('" + filename + "') failed.");
  }

};

// Combine two HMMs, using transition-matrix probability p1 for h1, and p2 for h2.
template <class T>
void combine( hmm_t<T> &H, const hmm_t<T> &h1, const hmm_t<T> &h2, const T p1 = 0.5, const T p2 = 0.5)
{
  if( h1.K != h2.K)
    throw std::runtime_error( "combine(): incompatible HMM state K's " + to_str(h1.K) + " and " + to_str(h2.K) + ".");
  H.K = h1.K;
  H.S = h1.S + h2.S;

  const size_t M  = h1.smps(0).dimensions();
  const size_t M2 = h2.smps(0).dimensions();
  if( M != M2)
    throw std::runtime_error( "combine(): incompatible HMM input sizes " + to_str(M) + " and " + to_str(M2) + ".");

  // Allocate memory.
  H.lPi.resize(H.S);
  H.lA .resize(H.S, H.S);
  H.smps.resize(H.S);
  for (int s=0; s<H.S; ++s) H.smps(s).init(M, H.K);

  // Concatenate h1.smps and h2.smps into H.smps, S-wise.
  // In h1, h2, and H, only S differs.  M and K are the same, by construction.
  for (int s=0; s<h1.S; ++s) H.smps(s       ).stuff(h1.smps(s));
  for (int s=0; s<h2.S; ++s) H.smps(s + h1.S).stuff(h2.smps(s));

  // Make transition matrix and initial probabilities:
  // concatenate h1.lPi and h2.lPi into H.lPi;
  // block-diagonal concatenate h1.lA and h2.lA into H.lA.
  for( size_t i=0; i < H.lA.size()/*SxS*/; ++i)
    H.lA[i] = -HUGE_VAL;
  for (int i=0; i < h1.S; ++i) {
    H.lPi(i) = h1.lPi(i);
    for (int j=0; j < h1.S; ++j) H.lA(i     ,j     ) = h1.lA(i,j);
  }
  for (int i=0; i < h2.S; ++i) {
    H.lPi(i+h1.S) = h2.lPi(i);
    for (int j=0; j < h2.S; ++j) H.lA(i+h1.S,j+h1.S) = h2.lA(i,j);
  }
  H.lA(h1.S-1       , h1.S) = log(p1);
  H.lA(h1.S-1 + h2.S,    0) = log(p2);

  // Normalize them.
  H.lPi.normalize_log();

  for (int s=0; s<h1.S; ++s) H.lA(h1.S-1     , s     ) += log1p(-p1);
  for (int s=0; s<h2.S; ++s) H.lA(h1.S-1+h2.S, s+h1.S) += log1p(-p2);
}

#endif
