// hmm.h -- Hidden Markov Model class with Gaussian Mixture Model states
// language: C++
// author  : Paris Smaragdis

#ifndef __HMM_H__
#define __HMM_H__

#include <fstream>
#include <algorithm>
#include <iostream>

#include "gmm.h"
#include "array.h"
#include "logadd.h"

// Hidden Markov Model class

template <class T>
class hmm_t {
public:
  int S; // States
  int K; // Gaussians per state
  array<T> lPi, lA; // Model parameters

  array<gmm_t<T> > gmms; // Gaussian mixture data
private:
  array<T> la, lb, xi, txi; // Parameters for Baum-Welch iterations

public:
  // Constructor
  hmm_t( int s = 0, int k = 1) : S(s), K(k) {}

  // Learn data
  void train( const array<T> &x, int iters = 100, const hmm_t<T> &H = hmm_t<T>())
  {
    // Input dims of array are reversed
    const int M = x.n;
    const int N = x.m;
    std::cout << "HMM training " << M << " x " << N << std::endl;

    // Setup state model parameters
    gmms.resize(S);
    for( int s=0; s<S; ++s){
      // M and N come from arg x.  K comes from constructor, or from load().
      gmms(s).ldt.resize(K);
      gmms(s).c.resize(K);
      gmms(s).m.resize(M,K);
      gmms(s).is.resize(M,K);
    }
    lPi.resize( S);
    lA.resize( S, S);

    if( H.S == 0){
      // Initial values of Gaussians
      for( int s=0; s<S; ++s){
	for( int k=0; k<K; ++k){
	  gmms(s).ldt(k) = 0;
	  gmms(s).c(k) = 1./K;
	  for( int i=0; i<M; ++i){
	    gmms(s).m(i,k) = T(rand())/RAND_MAX - 0.5;
	    const T init = 0.1;
	    gmms(s).is(i,k) = init;
	    gmms(s).ldt(k) += log(init);
	  }
	}
      }
      // Initial values of initial and transition probabilities
      for( int s=0; s<S; ++s)
	lPi(s) = log( 1./S);
      for( int i=0; i<S; ++i)
	for( int j=0; j<S; ++j)
	  lA(i,j) = log( T(rand())/RAND_MAX);
      for( int i=0; i<S; ++i){
	T ls = log(0.0);
	for( int j=0; j<S; ++j)
	  logadd( ls, lA(i,j));
	for( int j=0; j<S; ++j)
	  lA(i,j) -= ls;
      }
    }else{
      // Copy values of Gaussians
      gmms = H.gmms;

      // Copy values of initial and transition probabilities
      for( int i=0; i<S; ++i){
	lPi(i) = H.lPi(i);
	for( int j=0; j<S; ++j)
	  lA(i,j) = H.lA(i,j);
      }
    }

    // Allocate buffers
    array<T> q( N, K, S);
    array<T> g( N, K, S);
    array<T> lp( S, N);
    array<T> lk( iters);
    la .resize( S, N);
    lb .resize( S, N);
    xi .resize( S, S);
    txi.resize( S, S);

    // Iterate expectation-maximization
    for( int it=0; it<iters; ++it){

      // *** Expectation step ***

      // Get likelihoods from each gaussian from each state
#pragma omp parallel for
      for( int s=0; s<S; ++s){
	for( int k=0; k<K; ++k){
	  const T gc = log( gmms(s).c(k)) + 0.5*gmms(s).ldt(k) - 0.5*M*log(2*M_PI);
	  for( int j=0; j<N; ++j){
	    T qt = 0;
	    for (int i=0; i<M; ++i) {
	      qt += gmms(s).is(i,k) * sq(x(j,i) - gmms(s).m(i,k));
	    }
	    q(j,k,s) = gc - 0.5*qt;
	  }
	}
      }

      // Compute overall state likelihoods
      for( int s=0; s<S; ++s)
	for( int j=0; j<N; ++j){
	  T tp = log( 0.0);
	  for( int k=0; k<K; ++k)
	  logadd( tp, q(j,k,s));
	  lp(s,j) = tp;
	}

      // Get alphas
      for( int i=0; i<S; ++i)
	la(i,0) = lPi(i) + lp(i,0);
      for( int t=0; t<N-1; ++t)
	for( int j=0; j<S; ++j){
	  T ls = log( 0.0);
	  for( int i=0; i<S; ++i)
	    logadd( ls, la(i,t) + lA(i,j));
	  la(j,t+1) = ls + lp(j,t+1);
	}

      // Get betas
      for( int i=0; i<S; ++i)
	lb(i,N-1) = log( 0.0);
      lb(S-1,N-1) = 0;
      for( int t = N-2; t >= 0; --t)
	for( int i=0; i<S; ++i){
	  T ls = log( 0.0);
	  for( int j=0; j<S; ++j)
	    logadd( ls, lb(j,t+1) + lA(i,j) + lp(j,t+1));
	  lb(i,t) = ls;
	}

      // Get Xi
      for( int i=0; i<S*S; ++i)
	xi[i] = log(0.0);
      for( int t=0; t<N-1; ++t){
	T ls = log(0.0);
	for( int i=0; i<S; ++i)
	  for( int j=0; j<S; ++j){
	    txi(i,j) = lA(i,j) + la(i,t) + lp(j,t+1) + lb(j,t+1);
	    logadd( ls, txi(i,j));
	  }
	for( int i=0; i<S*S; ++i)
	  logadd( xi[i], txi[i] - ls);
      }

      // Get gamma
      for( int j=0; j<N; ++j){
	T ls = log( 0.0);
	for( int s=0; s<S; ++s)
	  logadd( ls, la(s,j)+lb(s,j));
	for( int s=0; s<S; ++s){
	  T tg = la(s,j) + lb(s,j) - lp(s,j) - ls;
	  for( int k=0; k<K; ++k)
	    g(j,k,s) = tg + q(j,k,s);
	}
      }

      // Get overall likelihood
      lk(it) = log( 0.0);
      for( int i=0; i<S; ++i)
	logadd( lk(it), la(i,N-1));
      if( ((it+1)%5==0) || it == iters-1)
	std::cout << "HMM iteration " << it+1 << " of " << iters << ": likelihood " << lk(it) << std::endl;

      // Exit log domain
      for( int i=0; i<N*K*S; ++i)
	g[i] = exp(g[i]);


      // *** Maximization step ***

      for (int s=0; s<S; ++s) {
	// Initial probabilities
	T tp = log(0.0);
	for( int i=0; i<K; ++i)
	  logadd(tp, log(g(0,i,s)));
	lPi(s) = tp;

	// Transition matrix
	T ls = log(0.0);
	for( int j=0; j<S; ++j)
	  logadd(ls, xi(s,j));
	for( int j=0; j<S; ++j)
	  lA(s,j) = xi(s,j) - ls;

	gmm_t<T>& gmm = gmms(s);

	// Priors
	for( int k=0; k<K; ++k){
	  T tc = 0;
	  for( int i=0; i<N; ++i)
	    tc += g(i,k,s);
	  gmm.c(k) = tc;
	}

	gmm.c.normalize();
      }
      // ******* IS SCALING RIGHT?  I get c = [1 1 1 1 1 ...]

#undef use_gmm_maximize
#ifdef use_gmm_maximize
      // TODO: instead of updating Means and Covariances in here, call gmm_t<T>::maximize(g, x, dummy);
      array<int> dummy(g.n);
      for (size_t _=0; _<dummy.size(); ++_)
	dummy[_] = true;
#endif

      for( int s=0; s<S; ++s){
	gmm_t<T>& gmm = gmms(s);
#ifdef use_gmm_maximize
	gmm.maximize(g, x, dummy);
#else
	// Means
	for( int k=0; k<K; ++k){
	  //const T ps = g.sum(k); // same as sg?
	  for( int i=0; i<M; ++i){
	    T ms = 0;
	    T sg = 0;
	    for( int j=0; j<N; ++j){
	      const T t = g(j,k,s);
	      sg += t;
	      ms += t * x(j,i);
	    }
	    gmm.m(i,k) = ms / sg;
	  }
	}

	// Covariances
	for( int k=0; k<K; ++k){
	  gmm.ldt(k) = 0;
	  for( int i=0; i<M; ++i){
	    T tu = 0;
	    T sg = 0;
	    for( int j=0; j<N; ++j){
	      const T t = g(j,k,s);
	      tu += t * sq(x(j,i)-gmms(s).m(i,k));
	      sg += t;
	    }
	    const T init = sg / tu;
	    gmm.is(i,k) = init;
	    gmm.ldt(k) += log(init);
	  }
	}
#endif
      }
    }
  }

  // Classify using a known HMM
  void classify( const array<T> &x, array<int> &q, const array<T> &bias = array<T>(), const int ist = -1)
  {
    // Input dims of array are reversed
    const int M = x.n, N = x.m;

    // Get state probabilities
    array<T> lB( S, N);
    for( int i=0; i<S*N; ++i)
      lB[i] = log(0.0); // Init all lB(s,j)'s.

#pragma omp parallel for
    for( int s=0; s<S; ++s)
      for( int k=0; k<K; ++k){
	const T gc = log(gmms(s).c(k)) + 0.5*gmms(s).ldt(k) - 0.5*M*log(2*M_PI);
	for( int j=0; j<N; ++j){
	  T qt = 0;
	  for( int i=0; i<M; ++i)
	    qt += gmms(s).is(i,k) * sq(x(j,i) - gmms(s).m(i,k));
	  logadd( lB(s,j), gc - 0.5*qt);
	}
      }

    // Add the bias to first state
    if (!bias.empty()) {
      if (int(bias.size()) != S)
	throw std::runtime_error( "classify(): bias vector has mismatched size");
      for (int i=0; i<S; ++i)
	for (int j=0; j<N; ++j)
	  lB(i,j) += bias(i);
    }

    viterbi( lB, q, ist);
  }

private:
  T sq(const T x) const { return x*x; }

  // Viterbi decoding
  void viterbi( const array<T> &lB, array<int> &q, const int ist = -1)
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

  // Short-time Viterbi decoding
  void stviterbi( const array<T> &lB, array<int> &q)
  {
    const int M = lB.m;
    const int N = lB.n;
    q.resize( N);

    // Initial probability vector
    array<T> pp( S);
    for( int i=0; i<S; ++i)
      pp(i) = lPi(i);

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
  void iviterbi( const array<T> &lB, const int a, const int b, array<T> &lP, array<int> &q)
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

public:
  void save( const std::string& filename)
  {
    using namespace std;
    if (filename.empty())
      throw runtime_error( "hmm_t::save(\"\") failed.");
    ofstream f( filename.c_str(), ios::out | ios::binary);
    if (!f)
      throw runtime_error( "hmm_t::save('" + filename + "') failed.");

    f.write((char*)&S,       sizeof(int)); // number of states
    f.write((char*)&lPi[0],  S*sizeof(T)); // initial log probabilities
    f.write((char*)&lA[0], S*S*sizeof(T)); // log transition matrix
    f.write((char*)&K,       sizeof(int)); // number of gaussians

    const int M = gmms[0].m.size() / K;
    std::cout << ( "DEBUG hmm_t::save('" + filename + "'): number of dimensions M is " + to_str(M) + ".\n");
    f.write((char*)&M, sizeof(int)); // dimension
    for (int s=0; s<S; ++s) {
      f.write((char*)&gmms[s] .c[0],   K*sizeof(T)); // priors
      f.write((char*)&gmms[s] .m[0], M*K*sizeof(T)); // means
      f.write((char*)&gmms[s].is[0], M*K*sizeof(T)); // inverse variances
    }

    if (!f)
      throw runtime_error( "hmm_t::save('" + filename + "') failed.");
    cout << "Saved HMM file " << filename << ".\n";
  }

  void load( const std::string& filename)
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
      throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of states, " + to_str(S) + ".");

    // initial log probabilities
    lPi.resize( S);
    f.read( (char*)&lPi[0], S*sizeof( T));

    // log transition matrix
    lA.resize( S, S);
    f.read( (char*)&lA[0], S*S*sizeof( T));

    // number of gaussians
    f.read( (char*)&K, sizeof( int));
    if (K <= 0)
      throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of gaussians per state, " + to_str(K) + ".");
    //std::cout << ( "DEBUG hmm_t::load('" + filename + "'): number of gaussians per state K is " + to_str(K) + ".\n");

    // dimension
    int M;
    f.read( (char*)&M, sizeof( int));
    if (M <= 0)
      throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of dimensions " + to_str(M) + ".");
    //std::cout << ( "DEBUG hmm_t::load('" + filename + "'): number of dimensions M is " + to_str(M) + ".\n");

    gmms.resize(S);
    for (int s=0; s<S; ++s) {
      gmms(s).ldt.resize(  K);
      gmms(s)  .c.resize(  K);
      gmms(s)  .m.resize(M,K);
      gmms(s) .is.resize(M,K);
    }
    for (int s=0; s<S; ++s) {
      f.read((char*)&gmms[s] .c[0],   K*sizeof(T)); // priors
      f.read((char*)&gmms[s] .m[0], M*K*sizeof(T)); // means
      f.read((char*)&gmms[s].is[0], M*K*sizeof(T)); // inverse variances
    }
    if (!f)
      throw runtime_error( "hmm_t::load('" + filename + "') failed.");

    // compute determinants
    for( int s=0; s<S; ++s) {
      for( int k=0; k<K; ++k) {
	T t=0;
	for( int i=0; i<M; ++i)
	  t += log(gmms(s).is(i,k));
	gmms(s).ldt(k) = t;
      }
    }
  }

};

// Combine two HMMs.
template <class T>
void combine( hmm_t<T> &H, const hmm_t<T> &h1, const hmm_t<T> &h2, T p1 = 0.5, T p2 = 0.5)
{
  if( h1.K != h2.K)
    throw std::runtime_error( "combine(): incompatible HMM state K's " + to_str(h1.K) + " and " + to_str(h2.K) + ".");
  H.K = h1.K;
  H.S = h1.S + h2.S;

  // Allocate parameters.  Second .m is array h.m's first dimension.
  if( h1.gmms(0).m.m != h2.gmms(0).m.m)
    throw std::runtime_error( "combine(): incompatible HMM input sizes " + to_str(h1.gmms(0).m.m) + " and " + to_str(h2.gmms(0).m.m) + ".");
  const size_t M = h1.gmms(0).m.m;

  H.lPi.resize(         H.S);
  H.lA .resize(    H.S, H.S);

  H.gmms.resize(H.S);
  for( int s=0; s<H.S; ++s){
    H.gmms(s).ldt.resize(   H.K);
    H.gmms(s).c  .resize(   H.K);
    H.gmms(s).m  .resize(M, H.K);
    H.gmms(s).is .resize(M, H.K);
  }

  // Concatenate h1.gmms and h2.gmms into H.gmms.
  // M and K are the same for both src's, and dst.  Only S differs.
  for( int s=0; s<h1.S; ++s){
    const gmm_t<T>& src = h1.gmms[s];
          gmm_t<T>& dst =  H.gmms[s];
    // Arrays m and is are 2D not 1D, but M and K are the same for
    // src and dst, so don't bother copying elementwise in a nested loop.
    std::copy(src.c  .v, src.c  .v + src.c  .size(), dst.c  .v);
    std::copy(src.ldt.v, src.ldt.v + src.ldt.size(), dst.ldt.v);
    std::copy(src.m  .v, src.m  .v + src.m  .size(), dst.m  .v);
    std::copy(src.is .v, src.is .v + src.is .size(), dst.is .v);
  }
  for( int s=0; s<h2.S; ++s){
    const gmm_t<T>& src = h2.gmms[s];
          gmm_t<T>& dst =  H.gmms[s + h1.S];
    std::copy(src.c  .v, src.c  .v + src.c  .size(), dst.c  .v);
    std::copy(src.ldt.v, src.ldt.v + src.ldt.size(), dst.ldt.v);
    std::copy(src.m  .v, src.m  .v + src.m  .size(), dst.m  .v);
    std::copy(src.is .v, src.is .v + src.is .size(), dst.is .v);
  }

  // Make transition matrix and initial probabilities.
  for( size_t i=0; i < H.lA.size(); ++i)
    H.lA[i] = -HUGE_VAL;
  for( int i=0; i < h1.S; ++i){
    H.lPi(i) = h1.lPi(i);
    for( int j=0; j < h1.S; ++j)
      H.lA(i,j) = h1.lA(i,j);
  }
  for( int i=0; i < h2.S; ++i){
    H.lPi(h1.S+i) = h2.lPi(i);
    for( int j=0; j < h2.S; ++j)
      H.lA(h1.S+i,h1.S+j) = h2.lA(i,j);
  }
  H.lA(h1.S        - 1, h1.S) = log(p1);
  H.lA(h1.S + h2.S - 1,    0) = log(p2);

  // Normalize them.
  H.lPi.normalize_log();

  for (int s=0; s<h1.S; ++s) H.lA(h1.S-1     , s     ) += log1p(-p1);
  for (int s=0; s<h2.S; ++s) H.lA(h1.S-1+h2.S, s+h1.S) += log1p(-p2);
}

#endif
