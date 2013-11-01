// gmm.h -- Gaussian Mixture Model class
// language: C++
// author  : Paris Smaragdis

#ifndef __GMM_H__
#define __GMM_H__

#include <fstream>
#include <iostream>

#include <x86intrin.h> // All SIMD intrinsics.  Requires -march=native, -mfpmath, -msse, etc; AND -On for n>0.

#include "array.h"
#include "logadd.h"
#include "state.h"

// Gaussian mixture model class
class gmm_t: public state_t {
  int K;		// number of gaussians
  int M;		// dimensions of input data
  array<real_t> c;	// priors
  array<real_t> m;	// means
  array<real_t> is;	// inverse variances
  array<real_t> ldt;	// covariances aka determinants
  real_t dg;		// Diagonal load
  // Always, M == m.m == is.m.

public:
  // Constructor
  gmm_t(int k=0/*uninitialized*/, real_t d=__FLT_EPSILON__) : K(k), m(0), dg(d) {}

  bool uninitialized() const { return K <= 0; }
  int dimensions()     const { return M; }
  int numGaussians()   const { return K; }

  void init(const int mArg, const int kArg)
  {
    M = mArg;
    K = kArg;
    c  .resize(  K);
    m  .resize(M,K);
    is .resize(M,K);
    ldt.resize(  K);
  }
  void init(const int mArg, const int kArg, std::ifstream& f)
  {
    init(mArg, kArg);
    readContents(f);
  }
  void init2()
  {
    if (uninitialized())
      throw std::runtime_error( "gmm_t::init2(): uninitialized.");
    const real_t startingValue = 0.1;
    for (int k=0; k<K; ++k) {
      ldt(k) = 0;
      c(k) = 1.0/K;
      for (int i=0; i<M; ++i) {
	m(i,k) = real_t(rand())/RAND_MAX - 0.5;
	is(i,k) = startingValue;
	ldt(k) += log(startingValue); // TODO: just assign ldt(k) directly, M * log(startingValue).
      }
    }
  }

  void normalizePriors()    { c.normalize(); }
  void normalizeLogPriors() { c.normalize_log_also(); }

  // Learn data "x"
  void train( const array<real_t> &x, const int iters = 100, const gmm_t& G = gmm_t(), bool prior = false)
  {
    if (uninitialized())
      throw std::runtime_error( "gmm_t::train(): uninitialized.");
    // Remember sizes
    const int M = x.n;
    const int N = x.m;

    // Check input
    for( int i=0; i<M*N; ++i)
      if(isinf(x[i]) || isnan(x[i]))
	throw std::runtime_error( "gmm_t::train() got infinity or NaN.");

    // Setup
    init(M, K);

    array<real_t> p( N, K);
    array<real_t> lk( iters);

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
	    real_t vm = 0;
	    for( int j = 0 ; j < N ; j++)
	    vm += x(j,i);
	    vm /= N;
	    real_t v = dg;
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
      likelihoods(p, x);

      // Massage posterior into shape and compute likelihood
      lk(it) = 0;
      for( int j = 0 ; j < N ; j++){
	real_t mx = p(j,0);
	for( int i = 1 ; i < K ; i++)
	  mx = std::max( mx, p(j,i));
	for( int i = 0 ; i < K ; i++)
	  p(j,i) -= mx;
	lk(it) += mx;
      }
      for( int j = 0 ; j < N ; j++){
	real_t t = p(j,0);
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
	array<real_t> m2(m), is2(is), c2(c), ldt2(ldt);
	array<int> learn2(learn);
	// Like init(M, nK), but without updating K yet.
	c  .resize(   nK);
	m  .resize(M, nK);
	is .resize(M, nK);
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
  void maximize(array<real_t>& p, const array<real_t> &x, const array<int>& learn)
  {
    if (uninitialized())
      throw std::runtime_error( "gmm_t::maximize(): uninitialized.");
    const int M = x.n;
    const int N = x.m;
    if (K != int(learn.size()) || K != int(p.n) || N != int(p.m))
      throw std::runtime_error( "gmm_t::maximize(): Incompatible sizes " + to_str(K) + " == " + to_str(learn.size()) + " == " + to_str(p.n) + "; " + to_str(N) + " == " + to_str(p.m) + ".");

#pragma omp parallel for
    for(int k=0; k<K; ++k){
      // Weights
      const real_t ps = p.sum(k);
      c(k) = ps/N;

      if (learn(k)) {
	// Means
	for (int i=0; i<M; ++i){
	  real_t ms = 0;
	  for (int j=0; j<N; ++j)
	    ms += x(j,i) * p(j,k);
	  // ms == cblas_ddot( N, &x.v[x.m*i], 1, &p.v[p.m*k], 1);
	  m(i,k) = ms/ps;
	}

	// Variances
	ldt(k) = 0;
	for (int i=0; i<M; ++i){
	  real_t ss = dg;
	  for (int j=0; j<N; ++j)
	    ss += sq(x(j,i)-m(i,k)) * p(j,k);
	  const real_t init = ps/ss;
	  is(i,k) = init;
	  ldt(k) += log(init);
	}
      }
    }
  }

  real_t placeholderForFunctionName(const array<real_t>& x, const int j, const int k) const
  {
    if (uninitialized())
      throw std::runtime_error( "gmm_t::placeholderForFunctionName(): uninitialized.");

    const real_t gc = log(c(k)) + 0.5*ldt(k) - 0.5*M*log(2*M_PI); // bug: move this out of the j-loops that call this
    real_t qt = 0;
    for (int i=0; i<M; ++i)
      qt += is(i,k) * sq(x(j,i) - m(i,k));
    return gc - 0.5*qt;
  }

  real_t placeholderForFunctionName2(const array<real_t>& x, const int j, const int k, const int kArg, const int mArg)
  {
    if (K<=0 && kArg>0) {
      K = kArg;
      M = mArg;
    }
    return placeholderForFunctionName(x,j,k);
  }

private:
  real_t sq(const real_t x) const { return x*x; }

  // Evaluate log likelihoods of M*N data x into N*K-array p.
  // Comparing several GMMs' likelihoods()s is like an HMM's classify().
  // (Maybe later, another bool arg disables error checking, for during training.)
  void likelihoods(array<real_t> &p, const array<real_t> &x)
  {
    if (uninitialized())
      throw std::runtime_error( "gmm_t::likelihoods(): uninitialized.");
    const int M = x.n;
    const int N = x.m;
    if (M != int(m.m))
      throw std::runtime_error( "gmm_t::likelihoods(): incompatible sizes.");
    if (N != int(p.m) || K != int(p.n))
      throw std::runtime_error( "gmm_t::likelihoods(): result array expected " + to_str(N) +" x "+ to_str(K) + ", not " + to_str(p.m) +" x "+ to_str(p.n) + ".");

    for (int k=0; k<K; ++k) {
      for (int j=0; j<N; ++j) {
	p(j,k) = placeholderForFunctionName(x,j,k);
      }
    }
  }

  void computeDeterminants()
  {
    for (int k=0; k<K; ++k) {
      real_t t=0;
      for (int i=0; i<M; ++i)
	t += log(is(i,k));
      ldt(k) = t;
    }
  }

public:
  // Read part of the data.  Compute the rest.
  void readContents(std::ifstream& f)
  {
    c .read(f, K); // priors
    m .read(f, M, K); // means
    is.read(f, M, K); // inverse variances
    computeDeterminants();
  }

  // Write part of the data.
  void writeContents(std::ofstream& f) const
  {
    c .write(f); // priors
    m .write(f); // means
    is.write(f); // inverse variances
  }

  void save(const std::string& filename) const
  {
    using namespace std;
    if (uninitialized())
      throw runtime_error( "gmm_t::save(): uninitialized.");
    if (filename.empty())
      throw runtime_error( "gmm_t::save(\"\") failed.");
    ofstream f( filename.c_str(), ios::out | ios::binary);
    if (!f)
      throw runtime_error( "gmm_t::save('" + filename + "') failed.");

    f.write((char*)&K, sizeof(int)); // number of gaussians
    f.write((char*)&M, sizeof(int)); // dimension
    writeContents(f);
    if (!f)
      throw runtime_error( "gmm_t::save('" + filename + "') failed.");
    cout << "Saved GMM file " << filename << ".\n";
  }

  void load(const std::string& filename)
  {
    using namespace std;
    ifstream f( filename.c_str(), ios::in | ios::binary);

    // number of gaussians
    f.read( (char*)&K, sizeof( int));
    if (K <= 0)
      throw std::runtime_error( "gmm_t::load(): nonpositive number of gaussians.");

    // dimension
    int M;
    f.read( (char*)&M, sizeof( int));
    if (M <= 0)
      throw std::runtime_error( "gmm_t::load(): nonpositive dimension.");

    init(M,K);
    readContents(f);
  }

  void copyGuts(const gmm_t& src)
  {
    c = src.c;
    ldt = src.ldt;
    m = src.m;
    is = src.is;
  }
  void stuff(const gmm_t& src)
  {
    // Arrays m and is are 2D not 1D, but M and K are the same for
    // src and dst when called from hmm_t::combine(),
    // so don't bother copying elementwise in a nested loop.
    std::copy(src.c  .v, src.c  .v + src.c  .size(), c  .v);
    std::copy(src.m  .v, src.m  .v + src.m  .size(), m  .v);
    std::copy(src.is .v, src.is .v + src.is .size(), is .v);
    std::copy(src.ldt.v, src.ldt.v + src.ldt.size(), ldt.v);
  }
};

#endif
