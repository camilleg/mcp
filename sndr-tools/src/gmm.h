// gmm.h -- Gaussian Mixture Model class
// language: C++
// author  : Paris Smaragdis

#ifndef __GMM_H__
#define __GMM_H__

#include <fstream>
#include <algorithm>
#include <iostream>

#include <x86intrin.h> // All SIMD intrinsics, not just SSE3.
// Has effect only with -march=native, -mfpmath, -msse, etc; AND -On for n>0.

#include "array.h"
#include "logadd.h"

// Gaussian mixture model class
template <class T>
class gmm_t {
public:
  // (Don't hide these publics behind accessors until the code stabilizes.)
  int K;        // Gaussians
  array<T> ldt;	// covariances aka determinants
  array<T> c;	// priors
  array<T> m;	// means
  array<T> is;	// inverse variances
private:
  T dg;  // Diagonal load

public:
  // Constructor with default values
  gmm_t( int k = 0, T d = __FLT_EPSILON__) : K( k), dg( d) {}

  // Learn data "x"
  void train( const array<T> &x, int iters = 100, const gmm_t<T> &G = gmm_t<T>(), bool prior = false)
  {
    // Remember sizes
    const int M = x.n;
    const int N = x.m;

    // Check input
    for( int i=0; i<M*N; ++i)
      if(isinf(x[i]) || isnan(x[i]))
	throw std::runtime_error( "gmm_t::train() got infinity or NaN.");

    // Setup
    ldt.resize( K);
    c.resize( K);
    m.resize( M, K);
    is.resize( M, K);

    array<T> p( N, K);
    array<T> lk( iters);

#ifndef __NO_RAND_SEED
    // Randomize
    time_t tm; time( &tm);
    srand( tm);
    (void)rand();
#endif

    // Sort out the learning situation
    array<int> learn( K);
    if( G.K != 0 && prior){
      for( int k = 0 ; k < K ; ++k)
	learn(k) = k >= G.K;
    }else
      for( int k = 0 ; k < K ; ++k)
	learn(k) = 1;

    // Initial values
    for( int k = 0 ; k < K ; k++){
	if( learn( k)){
	  ldt(k) = 0;
	  c(k) = 1.;
	  for( int i = 0 ; i < M ; i++){
	    const int ri = (N-1)*double( rand())/RAND_MAX;
	    m(i,k) = x(ri,i);
	    //					m(i,k) = x(k,i);
	    //					m(i,k) = T( rand())/RAND_MAX - .5;
	    T vm = 0;
	    for( int j = 0 ; j < N ; j++)
	    vm += x(j,i);
	    vm /= N;
	    T v = dg;
	    for( int j = 0 ; j < N ; j++)
	    v += (x(j,i)-vm)*(x(j,i)-vm);
	    v = v / (N-1);
	    is(i,k) = 1./(v+dg);
	    //					is(i,k) = .1;1./1000;
	    ldt(k) += log( is(i,k));
	  }
	}else{
	  ldt(k) = G.ldt(k);
	  c(k) = G.c(k);
	  for( int i = 0 ; i < M ; i++){
	    m(i,k) = G.m(i,k);
	    is(i,k) = G.is(i,k);
	  }
	}
      }
      c.normalize();

    // Iterate expectation-maximization.
    for( int it = 0 ; it < iters ; it++){

      // *** Expectation step ***

#pragma omp parallel for
      for( int k = 0 ; k < K ; k++){
	T gc = log( c(k)) + 0.5*ldt(k) - 0.5*M*log(2*M_PI);
	for( int j = 0 ; j < N ; j++){
	  T qt = 0;
	  for( int i = 0 ; i < M ; i++)
	    qt += is(i,k) * (x(j,i) - m(i,k)) * (x(j,i) - m(i,k));
	  p(j,k) = gc - 0.5*qt;
	}
      }

      // Massage posterior into shape and compute likelihood
      lk(it) = 0;
      for( int j = 0 ; j < N ; j++){
	T mx = p(j,0);
	for( int i = 1 ; i < K ; i++)
	  mx = std::max( mx, p(j,i));
	for( int i = 0 ; i < K ; i++)
	  p(j,i) -= mx;
	lk(it) += mx;
      }
      for( int j = 0 ; j < N ; j++){
	T t = p(j,0);
	for( int i = 1 ; i < K ; i++)
	  logadd(t, p(j,i));
	for( int i = 0 ; i < K ; i++)
	  p(j,i) -= t;
	lk(it) += t;
      }
      if( ((it+1)%5==0) || it == iters-1)
	std::cout << "GMM iteration " << it+1 << " of " << iters << ": likelihood " << lk(it) << std::endl;

      // Exit log domain
      for (int i=0; i<K*N; ++i)
	p[i] = exp(p[i]);

      // *** Maximization step ***
      maximize(p, x, learn);

      // Remove blown-up states
      int nK = K;
      for (int k=0; k<K; ++k)
	if (isinf(ldt(k)) || isnan(ldt(k)))
	  --nK;
      if (nK != K){
	std::cout << "GMM iteration " << it << ", removing " << K-nK << " blown-up states." << std::endl;
	// Copy constructor is overkill, because the .v[] std::copy's will get clobbered in the k-loop.
	// But at least the copy constructor correctly sets dimensions of the xxx2 arrays.
	array<T> m2(m), is2(is), c2(c), ldt2(ldt);
	array<int> learn2(learn);
	m  .resize(M, nK);
	is .resize(M, nK);
	c  .resize(   nK);
	ldt.resize(   nK);
	p.resize(N, nK);
	for (int k=0, ck=0; k<K; ++k){
	  if (!isinf(ldt2(k)) && !isnan(ldt2(k))) {
	    // Copy nonblownup state from old arrays ("2", k, K) to new arrays (ck, nK).
	    learn(ck) = learn2(k);
	    ldt(ck) = ldt2(k);
	    c(ck) = c2(k);
	    for (int i=0; i<M; ++i) {
	      m (i,ck) = m2 (i,k);
	      is(i,ck) = is2(i,k);
	    }
	    ++ck;
	  }
	}
	K = nK;
      }
    }
  }

  // Maximization step of expectation-maximization.
  // Update arg p and members c, m, is, ldt.
  void maximize(array<T>& p, const array<T> &x, const array<int>& learn)
  {
    const int M = x.n;
    const int N = x.m;
    if (K != int(learn.size()) || K != int(p.n) || N != int(p.m))
      throw std::runtime_error( "gmm_t::maximize(): Incompatible sizes");
#pragma omp parallel for
    for(int k=0; k<K; ++k){
      // Weights
      const T ps = p.sum(k);
      c(k) = ps/N;

      if (learn(k)) {
	// Means
	for (int i=0; i<M; ++i){
	  T ms = 0;
	  for (int j=0; j<N; ++j)
	    ms += x(j,i) * p(j,k);
	  // ms == cblas_ddot( N, &x.v[x.m*i], 1, &p.v[p.m*k], 1);
	  m(i,k) = ms/ps;
	}

	// Variances
	ldt(k) = 0;
	for (int i=0; i<M; ++i){
	  T ss = dg;
	  for (int j=0; j<N; ++j)
	    ss += (x(j,i)-m(i,k)) * (x(j,i)-m(i,k)) * p(j,k);
	  is(i,k) = ps/ss;
	  ldt(k) += log( is(i,k));
	}
      }
    }
  }

private:
  // Evaluate log likelihoods of data x into p.
  // Comparing several GMMs' likelihoods()s is like an HMM's classify().
  void likelihoods( const array<T> &x, array<T> &p)
  {
    // Check sizes and allocate output
    const int M = x.n;
    const int N = x.m;
    if( M != m.m)
      throw std::runtime_error( "gmm_t::likelihoods(): Incompatible sizes");
    p.resize( N);
    for( int i = 0 ; i < N ; i++)
      p(i) = log( 0.);

    //#pragma omp parallel for
    for( int k = 0 ; k < K ; k++){
      const T gc = log( c(k)) + 0.5*ldt(k) - 0.5*M*log(2*M_PI);
      for( int j = 0 ; j < N ; j++){
	T qt = 0;
	for( int i = 0 ; i < M ; i++)
	  qt += is(i,k) * (x(j,i) - m(i,k)) * (x(j,i) - m(i,k));
	logadd( p(j), gc - 0.5*qt);
      }
    }
  }

public:
  void save( const std::string& filename)
  {
    using namespace std;
    if (filename.empty())
      throw runtime_error( "gmm_t::save(\"\") failed.");
    ofstream f( filename.c_str(), ios::out | ios::binary);
    if (!f)
      throw runtime_error( "gmm_t::save('" + filename + "') failed.");

    f.write((char*)&K,          sizeof(int)); // number of gaussians
    f.write((char*)&m.m,        sizeof(int)); // dimension
    f.write((char*)&c[0],       K*sizeof(T)); // priors
    f.write((char*)&m[0],   m.m*K*sizeof(T)); // means
    f.write((char*)&is[0], is.m*K*sizeof(T)); // inverse variances
    if (!f)
      throw runtime_error( "gmm_t::save('" + filename + "') failed.");
    cout << "Saved GMM file " << filename << ".\n";
  }

  void load( const std::string& filename)
  {
    using namespace std;
    ifstream f( filename.c_str(), ios::in | ios::binary);

    // number of gaussians
    f.read( (char*)&K, sizeof( int));
    if( K <= 0)
      throw std::runtime_error( "gmm_t::load(): nonpositive number of gaussians.");

    // dimension
    int M;
    f.read( (char*)&M, sizeof( int));
    if( M <= 0)
      throw std::runtime_error( "gmm_t::load(): nonpositive dimension.");

    // priors
    c.resize(K);
    f.read((char*)&c[0], K*sizeof(T));

    // means
    m.resize(M, K);
    f.read((char*)&m[0], M*K*sizeof(T));

    // inverse variances
    is.resize(M, K);
    f.read((char*)&is[0], M*K*sizeof(T));

    // compute determinants
    ldt.resize( K);
    for( int k = 0 ; k < K ; k++){
      ldt(k) = 0;
      for( int i = 0 ; i < M ; i++)
	ldt(k) += log( is(i,k));
    }
  }

};

#endif
