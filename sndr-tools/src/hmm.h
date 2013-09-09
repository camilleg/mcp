// hmm.h -- Hidden Markov Model class with Gaussian Mixture Model states
// language: C++
// author  : Paris Smaragdis

#ifndef __HMM_H__
#define __HMM_H__

#include <cmath>
#include <fstream>
#include <algorithm>
#include <iostream>

#include "array.h"

// Hidden Markov Model class

template <class T>
class hmm_t {
public:
	int S; // States
	int K; // Gaussians per state
	array<T> lPi, lA; // Model parameters
	array<T> ldt, c, m, is; // Gaussian mixture data
	array<T> la, lb, xi, txi; // Various parameters for Baum-Welch iterations

	// Constructor
	hmm_t( int s = 0, int k = 1) : S(s), K(k) {}

	// Log-add y to x.  For log-sum-exp idiom, to avoid underflow and overflow.
	inline void logadd( T& x, const T& y)
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
	    throw std::runtime_error( "hmm_t::logadd(" + to_str(x) + ", " + to_str(y) + ")");

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

	// Learn data
	void train( const array<T> &x, int iters = 100, const hmm_t<T> &H = hmm_t<T>())
	{
		// Input dims of array are reversed
		const int M = x.n, N = x.m;
		std::cout << "HMM training " << M << " x " << N << std::endl;

		// Setup state model parameters
		ldt.resize( K, S);
		c.resize( K, S);
		m.resize( M, K, S);
		is.resize( M, K, S);
		lPi.resize( S);
		lA.resize( S, S);

		// Initial values
		if( H.S == 0){

			// Initial values for Gaussians
			for( int s = 0; s < S; ++s){
				for( int k = 0; k < K; ++k){
					ldt(k,s) = 0;
					c  (k,s) = 1./K;
					for( int i = 0; i < M; ++i){
						m (i,k,s) = T( rand())/RAND_MAX - 0.5;
						is(i,k,s) = 0.1;
						ldt(k,s) += log( is(i,k,s));
					}
				}
			}

			// Initial values for initial and transition probabilities
			for( int s = 0; s < S; ++s)
				lPi(s) = log( 1./S);
			for( int i = 0; i < S; ++i)
				for( int j = 0; j < S; ++j)
					lA(i,j) = log( T(rand())/RAND_MAX);
			for( int i = 0; i < S; ++i){
				T ls = log(0.0);
				for( int j = 0; j < S; ++j)
					logadd( ls, lA(i,j));
				for( int j = 0; j < S; ++j)
					lA(i,j) -= ls;
			}
		}else{
			
			// Copy values of Gaussians
			for( int s = 0; s < S; ++s){
				for( int k = 0; k < K; ++k){
					ldt(k,s) = H.ldt(k,s);
					c(k,s) = H.c(k,s);
					for( int i = 0; i < M; ++i){
						m (i,k,s) = H.m (i,k,s);
						is(i,k,s) = H.is(i,k,s);
					}
				}
			}

			// Copy values of initial and transition probabilities
			for( int i = 0; i < S; ++i){
				lPi(i) = H.lPi(i);
				for( int j = 0; j < S; ++j)
					lA(i,j) = H.lA(i,j);
			}
		}

		// Allocate buffers
		array<T> lp( S, N);
		array<T> q( N, K, S);
		array<T> g( N, K, S);
		la.resize( S, N);
		lb.resize( S, N);
		xi.resize( S, S);
		txi.resize( S, S);
		array<T> lk( iters);

		// Iterate expectation-maximization
		for( int it = 0; it < iters; ++it){
		
			// *** Expectation step ***

			// Get likelihoods from each gaussian from each state
			for( int s = 0 ; s < S ; s++){
				for( int k = 0 ; k < K ; k++){
					T gc = log( c(k,s)) + 0.5*ldt(k,s) - 0.5*M*log(2*M_PI);
					for( int j = 0 ; j < N ; j++){
						T qt = 0;
						for( int i = 0 ; i < M ; i++)
							qt += is(i,k,s) * (x(j,i) - m(i,k,s)) * (x(j,i) - m(i,k,s));
						q(j,k,s) = gc - 0.5*qt;
					}
				}
			}

			// Get overall state likelihoods
			for( int s = 0 ; s < S ; s++)
				for( int j = 0 ; j < N ; j++){
					T tp = log( 0.0);
					for( int k = 0 ; k < K ; k++)
						logadd( tp, q(j,k,s));
					lp(s,j) = tp;
				}

			// Get alphas
			for( int i = 0 ; i < S ; i++)
				la(i,0) = lPi(i) + lp(i,0);
			for( int t = 0 ; t < N-1 ; t++)
				for( int j = 0 ; j < S ; j++){
					T ls = log( 0.0);
					for( int i = 0 ; i < S ; i++)
						logadd( ls, la(i,t) + lA(i,j));
					la(j,t+1) = ls + lp(j,t+1);
				}

			// Get betas
			for( int i = 0 ; i < S ; i++)
				lb(i,N-1) = log( 0.0);
			lb(S-1,N-1) = 0;
			for( int t = N-2 ; t > -1 ; t--)
				for( int i = 0 ; i < S ; i++){
					T ls = log( 0.0);
					for( int j = 0 ; j < S ; j++)
						logadd( ls, lb(j,t+1) + lA(i,j) + lp(j,t+1));
					lb(i,t) = ls;
				}

			// Get Xi
			for( int i = 0 ; i < S*S ; i++)
				xi(i) = log( 0.0);
			for( int t = 0 ; t < N-1 ; t++){
				T ls = log( 0.0);
				for( int i = 0 ; i < S ; i++)
					for( int j = 0 ; j < S ; j++){
						txi(i,j) = lA(i,j) + la(i,t) + lp(j,t+1) + lb(j,t+1);
						logadd( ls, txi(i,j));
					}
				for( int i = 0 ; i < S*S ; i++)
					logadd( xi(i), txi(i) - ls);
			}

			// Get gamma
			for( int j = 0 ; j < N ; j++){
				T ls = log( 0.0);
				for( int s = 0 ; s < S ; s++)
					logadd( ls, la(s,j)+lb(s,j));
				for( int s = 0 ; s < S ; s++){
					T tg = la(s,j) + lb(s,j) - lp(s,j) - ls;
					for( int k = 0 ; k < K ; k++)
						g(j,k,s) = tg + q(j,k,s);
				}
			}

			// Get overall likelihood
			lk(it) = log( 0.0);
			for( int i = 0 ; i < S ; i++) {
				logadd( lk(it), la(i,N-1));
			}
			if( 1){ //( !((it+1)%25) | it == iters-1){
				std::cout << "HMM iteration " << it+1 << " of " << iters << ": likelihood " << lk(it) << std::endl;
			}

			// Get out of log domain
			for( int i = 0 ; i < N*K*S ; i++)
				g(i) = exp( g(i));

			// *** Maximization step ***

			// Initial probabilities
			for( int s = 0 ; s < S ; s++){
				T tp = log( 0.0);
				for( int i = 0 ; i < K ; i++)
					logadd( tp, log( g(0,i,s)));
				lPi(s) = tp;
			}

			// Transition matrix
			for( int i = 0 ; i < S ; i++){
				T ls = log( 0.0);
				for( int j = 0 ; j < S ; j++)
					logadd( ls, xi(i,j));
				for( int j = 0 ; j < S ; j++)
					lA(i,j) = xi(i,j) - ls;
			}

			// Priors
			for( int s = 0 ; s < S ; s++){
				for( int k = 0 ; k < K ; k++){
					T tc = 0;
					for( int i = 0 ; i < N ; i++)
						tc += g(i,k,s);
					c(k,s) = tc;
				}
			}
			for( int s = 0 ; s < S ; s++){
				T sc = 0;
				for( int k = 0 ; k < K ; k++)
					sc += c(k,s);
				for( int k = 0 ; k < K ; k++)
					c(k,s) /= sc;
			}
			// ******* IS SCALING RIGHT?  I get c = [1 1 1 1 1 ...]

			// Means
			for( int s = 0 ; s < S ; s++){
				for( int k = 0 ; k < K ; k++){
					for( int i = 0 ; i < M ; i++){
						T ms = 0, sg = 0;
						for( int j = 0 ; j < N ; j++){
							ms += x(j,i) * g(j,k,s);
							sg += g(j,k,s);
						}
						m(i,k,s) = ms / sg;
					}
				}
			}

			// Covariances
			for( int s = 0 ; s < S ; s++){
				for( int k = 0 ; k < K ; k++){
					ldt(k,s) = 0;
					for( int i = 0 ; i < M ; i++){
						T tu = 0, sg = 0;
						for( int j = 0 ; j < N ; j++){
							tu += (x(j,i)-m(i,k,s))*(x(j,i)-m(i,k,s))*g(j,k,s);
							sg += g(j,k,s);
						}
						is(i,k,s) = sg / tu;
						ldt(k,s) += log( is(i,k,s));
					}
				}
			}
		}
	}

	// Classify using a known HMM
	void classify( const array<T> &x, array<int> &q, array<T> &bias = array<T>(), int ist = -1)
	{
		// Input dims of array are reversed
		const int M = x.n, N = x.m;

		// Get state probabilities
		array<T> lB( S, N);
		for( int i = 0 ; i < S*N ; ++i)
			lB(i) = -log(0.0); // Init all lB(s,j)'s.

#pragma omp parallel for
		for( int s = 0 ; s < S ; s++)
			for( int k = 0 ; k < K ; k++){
				const T gc = log( c(k,s)) + 0.5*ldt(k,s) - 0.5*M*log(2*M_PI);
				for( int j = 0 ; j < N ; j++){
					T qt = 0;
					for( int i = 0 ; i < M ; i++)
						qt += is(i,k,s) * (x(j,i) - m(i,k,s)) * (x(j,i) - m(i,k,s));
					logadd( lB(s,j), gc - 0.5*qt);
				}
			}

		// Add the bias to first state
		if( bias.size()){
			if( int(bias.size()) != S)
				throw std::runtime_error( "clasify(): Bias vector has wrong size");
			for( int i = 0 ; i < S ; i++)
				for( int j = 0 ; j < N ; j++)
					lB(i,j) += bias(i);
		}

		// Run viterbi
		viterbi( lB, q, ist);
	}

	// Viterbi decoding
	void viterbi( const array<T> &lB, array<int> &q, const int ist = -1)
	{
		const int N = lB.n;

		// Temp space
		array<T> d( S, 2);
		array<int> p( S, N);

		// Override initial state if one is provided
		array<double> nlPi(S);
		if( ist != -1){
			for( int i = 0 ; i < S ; i++)
				nlPi(i) = log( 0.0);
			nlPi(ist) = 0.0;
		}else
			for( int i = 0 ; i < S ; i++)
				nlPi(i) = lPi(i);

		// Initialize
		int di = 0;
		for( int i = 0 ; i < S ; i++){
			d(i,di) = nlPi(i) + lB(i,0);
			p(i,0) = 0;
		}

		// Propagate
		for( int t = 1 ; t < N ; t++){
			for( int i = 0 ; i < S ; i++){

				// Find max likelihood and its index on temp variables
				T dt = -HUGE_VAL;
				int pt = 0;
				for( int j = 0 ; j < S ; j++){
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

		// Allocate output
		q.resize( N);

		// Terminate
		T l = d(0,di);
		q(N-1) = 0;
		for( int i = 1 ; i < S ; i++)
			if( d(i,di) > l){
				l = d(i,di);
				q(N-1) = i;
			}

		// Backtrack
		for( int t = N-2 ; t >= 0 ; t--)
			q(t) = p(q(t+1),t+1);
	}

	// Short-time Viterbi decoding
	void stviterbi( const array<T> &lB, array<int> &q)
	{
		const int M = lB.m;
		const int N = lB.n;

		// Allocate output
		q.resize( N);

		// Initial probability vector
		array<T> pp( S);
		for( int i = 0 ; i < S ; i++)
			pp(i) = lPi(i);

		int a = 0;
		for( int b=1; b<N; ){

			// Get state paths for a short segment
			array<int> s;
			iviterbi( lB, a, b, pp, s);

			// Find fusion point (where all possible solutions converge to the same state)
			int t = 0;
			for( int j = 0 ; j < s.n ; j++){
				bool e = true;
				for( int i = 0 ; i < S-1 ; i++)
					e &= s(i,j) == s(i+1,j);

				if( !e){
					t = j;
					break;
				}
			}

			if( t > 0){
				// Keep advance
				for( int i = 0 ; i < t ; i++)
					q(a+i) = s(0,i);
				for( int i = 0 ; i < S ; i++)
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

		// Temp space
		array<T> d( S, 2);
		array<int> p( S, N);
		
		// Initialize
		int di = 0;
		for( int i = 0 ; i < S ; i++){
			d(i,di) = lP(i) + lB(i,a);
			p(i,0) = 0;
		}

		// Propagate
		for( int t = 1 ; t < N ; t++){
			for( int i = 0 ; i < S ; i++){
				
				// Find max likelihood and its index on temp variables
				T dt = -HUGE_VAL;
				int pt = 0;
				for( int j = 0 ; j < S ; j++){
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

		// Allocate output
		q.resize( S, N);

		// Try all possible terminations
		for( int i = 0 ; i < S ; i++){
			q(i,N-1) = i;
			
			// Backtrack
			for( int t = N-2 ; t >= 0 ; t--)
				q(i,t) = p(q(i,t+1),t+1);
		}
	}

	// Save model
	void save( const std::string& filename)
	{
		using namespace std;
		if (filename.empty()) {
		  throw runtime_error( "hmm_t::save(\"\") failed.");
		}
		ofstream f( filename.c_str(), ios::out | ios::binary);
		if (!f) {
		  throw runtime_error( "hmm_t::save('" + filename + "') failed.");
		}

		// number of states
		f.write( (char*)&S, sizeof( int));

		// initial log probabilities
		f.write( (char*)&lPi(0), S*sizeof( T));

		// log transition matrix
		f.write( (char*)&lA(0), S*S*sizeof( T));

		// number of gaussians
		f.write( (char*)&K, sizeof( int));

		// dimension
		const int M = m.size()/(K*S);
		f.write( (char*)&M, sizeof( int));

		// priors
		f.write( (char*)&c(0), K*S*sizeof( T));

		// means
		f.write( (char*)&m(0), M*K*S*sizeof( T));

		// inverse variances
		f.write( (char*)&is(0), M*K*S*sizeof( T));

		if (!f) {
		  throw runtime_error( "hmm_t::save('" + filename + "') failed.");
		}
		cout << "Saved HMM file " << filename << ".\n";
	}

	// Load model
	void load( const std::string& filename)
	{
		using namespace std;
		if (filename.empty()) {
		  throw runtime_error( "hmm_t::load(\"\") failed.");
		}
		ifstream f( filename.c_str(), ios::in | ios::binary);
		if (!f) {
		  throw runtime_error( "hmm_t::load('" + filename + "') failed.");
		}
		
		// number of states
		f.read( (char*)&S, sizeof( int));
		if (S <= 0) {
		  throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of states, " + to_str(S) + ".");
		}

		// initial log probabilities
		lPi.resize( S);
		f.read( (char*)&lPi(0), S*sizeof( T));

		// log transition matrix
		lA.resize( S, S);
		f.read( (char*)&lA(0), S*S*sizeof( T));

		// number of gaussians
		f.read( (char*)&K, sizeof( int));
		if (K <= 0) {
		  throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of gaussians per state, " + to_str(K) + ".");
		}

		// dimension
		int M;
		f.read( (char*)&M, sizeof( int));
		if (M <= 0) {
		  throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive number of dimensions, from scalars " + to_str(m.size()) + ", " + to_str(K) + ", " + to_str(S) + ".");
		}

		// priors
		c.resize( K, S);
		f.read( (char*)&c(0), K*S*sizeof( T));

		// means
		m.resize( M, K, S);
		f.read( (char*)&m(0), M*K*S*sizeof( T));

		// inverse variances
		is.resize( M, K, S);
		f.read( (char*)&is(0), M*K*S*sizeof( T));

		if (!f) {
		  throw runtime_error( "hmm_t::load('" + filename + "') failed.");
		}
		
		// number of states
		f.read( (char*)&S, sizeof( int));
		if (S <= 0 || S > 1e9) {
		  throw runtime_error( "hmm_t::load('" + filename + "'): nonpositive or suspicious number of states, " + to_str(S) + ".");
		}

		// Compute the determinants.
		ldt.resize( K, S);
		for( int s = 0 ; s < S ; ++s)
			for( int k = 0 ; k < K ; ++k){
				ldt(k,s) = 0;
				for( int i = 0 ; i < M ; ++i)
					ldt(k,s) += log( is(i,k,s));
			}
	}

};

// Combine two HMMs in a loop
template <class T>
void combine( hmm_t<T> &H, const hmm_t<T> &h1, const hmm_t<T> &h2, T p1 = 0.5, T p2 = 0.5)
{
	H.S = h1.S + h2.S;
	H.K = h1.K;
	if( H.K != h2.K)
		throw std::runtime_error( "combine(): HMM state K's are incompatible");

	// Allocate parameters.  Second .m is array h.m's first dimension.
	const size_t M = h1.m.m;
	if( M != h2.m.m)
		throw std::runtime_error( "combine(): Input sizes are incompatible");
	H.lPi.resize(         H.S);
	H.lA .resize(    H.S, H.S);
	H.c  .resize(    H.K, H.S);
	H.ldt.resize(    H.K, H.S);
	H.m  .resize( M, H.K, H.S);
	H.is .resize( M, H.K, H.S);

	// Copy data
	for( int s = 0; s < h1.S; ++s)
		for( int k = 0; k < h1.K; ++k){
			H.c  (k,s) = h1.c  (k,s);
			H.ldt(k,s) = h1.ldt(k,s);
			for( size_t i = 0; i < M; ++i){
				H.m (i,k,s) = h1.m (i,k,s);
				H.is(i,k,s) = h1.is(i,k,s);
			}
		}
	for( int s = 0; s < h2.S; ++s)
		for( int k = 0; k < h2.K; ++k){
			H.c  (k,h1.S+s) = h2.c  (k,s);
			H.ldt(k,h1.S+s) = h2.ldt(k,s);
			for( size_t i = 0; i < M; ++i){
				H.m (i,k,h1.S+s) = h2.m (i,k,s);
				H.is(i,k,h1.S+s) = h2.is(i,k,s);
			}
		}
	
	// Make transition matrix and initial probabilities
	for( size_t i = 0; i < H.lA.size(); ++i)
		H.lA(i) = -HUGE_VAL;
	for( int i = 0; i < h1.S; ++i){
		H.lPi(i) = h1.lPi(i);
		for( int j = 0; j < h1.S; ++j)
			H.lA(i,j) = h1.lA(i,j);
	}
	for( int i = 0; i < h2.S; ++i){
		H.lPi(h1.S+i) = h2.lPi(i);
		for( int j = 0; j < h2.S; ++j)
			H.lA(h1.S+i,h1.S+j) = h2.lA(i,j);
	}
	H.lA(h1.S-1,      h1.S) = log( p1);
	H.lA(h1.S+h2.S-1, 0)    = log( p2);

	// Normalize them
	T ps = 0;
	for( int i = 0; i < H.S; ++i)
		ps += exp( H.lPi(i));
	for( int i = 0; i < H.S; ++i)
		H.lPi(i) = log( exp( H.lPi(i))/ps);

	for( int i = 0; i < h1.S; ++i)
		H.lA(h1.S-1,i) += log( (1.-p1));
	for( int i = 0; i < h2.S; ++i)
		H.lA(h1.S+h2.S-1,h1.S+i) += log( (1.-p2));
}

#endif
