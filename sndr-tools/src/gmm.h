// gmm.h -- Gaussian Mixture Model class
// language: C++
// author  : Paris Smaragdis

#ifndef __GMM_H__
#define __GMM_H__

#include <cmath>
#include <fstream>
#include <algorithm>
#include <iostream>

#include <x86intrin.h> // All SIMD intrinsics, not just SSE3.
// Has effect only with -march=native, -mfpmath, -msse, etc; AND -On for n>0.

#include "array.h"

// Gaussian mixture model class
template <class T>
class gmm_t {
public:
	int K; // Gaussians
private:
	T dg;  // Diagonal load
public:
	array<T> ldt;	// ?
	array<T> c;	// priors
	array<T> m;	// means
	array<T> is;	// inverse variances

	// Constructor with default values
	gmm_t( int k = 0, T d = __FLT_EPSILON__) : K( k), dg( d) {}
private:

	// log sum
	inline T lsum( T x, T y)
	{
		return x == y ?
			x+log(2.) :
			std::max( x, y) + log1p( exp( -fabs( x-y)));
	}

public:
	// Learn data
	void train( const array<T> &x, int iters = 100, const gmm_t<T> &G = gmm_t<T>(), bool prior = false)
	{
		// Remember sizes
		const int M = x.n;
		const int N = x.m;

		// Check input
		for( int i = 0 ; i < M*N ; ++i)
			if( isinf( x(i)) || isnan( x(i)))
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

		//  Normalize priors
		{
			T cs = 0;
			for( int k = 0 ; k < K ; k++)
				cs += c(k);
			for( int k = 0 ; k < K ; k++)
				c(k) /= cs;
		}

		// Start iterating
		for( int it = 0 ; it < iters ; it++){

			//
			// Expectation step
			//

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
				T mx = p(j);
				for( int i = 1 ; i < K ; i++)
					mx = std::max( mx, p(j,i));
				for( int i = 0 ; i < K ; i++)
					p(j,i) -= mx;
				lk(it) += mx;
			}
			for( int j = 0 ; j < N ; j++){
				T t = p(j);
				for( int i = 1 ; i < K ; i++)
					t = lsum( t, p(j,i));
				for( int i = 0 ; i < K ; i++)
					p(j,i) -= t;
				lk(it) += t;
			}
			if( !(it%25) || it == iters-1){
			  std::cout << "GMM Iteration " << it+1 << " of " << iters << ": likelihood " << lk(it) << std::endl;
			}

			// Get out of log domain
			for( int i = 0 ; i < K*N ; i++)
				p(i) = exp( p(i));


			//
			// Maximization step
			//

#pragma omp parallel for
			for( int k = 0 ; k < K ; k++){
				T ps = 0;
				for( int i = 0 ; i < N ; i++)
					ps += p(i,k);

				// Weights
				c(k) = ps/N;

				if( learn(k)) {
					// Means
					for( int i = 0 ; i < M ; i++){
						T ms = 0;
						for( int j = 0 ; j < N ; j++)
							ms += x(j,i) * p(j,k);
						m(i,k) = ms/ps;
					//	m(i,k) = cblas_ddot( N, &x.v[x.m*i], 1, &p.v[p.m*k], 1)/ps;
					}

					// Variances
					ldt(k) = 0;
					for( int i = 0 ; i < M ; i++){
						T ss = dg;
						for( int j = 0 ; j < N ; j++)
							ss += (x(j,i)-m(i,k))*(x(j,i)-m(i,k))*p(j,k);
						is(i,k) = ps/ss;
						ldt(k) += log( is(i,k));
					}
				}
			}

			int nK = K;
			for( int k=0; k<K; ++k)
				if( isinf( ldt(k)) || isnan( ldt(k)) )
					--nK;
			if( nK != K){
				std::cout << "GMM iteration " << it << ", " << K-nK << " removing blown-up states" << std::endl;
				array<T> m2( m), is2( is), c2( c), ldt2( ldt);
				array<int> learn2( learn);
				m.resize( M, nK);
				is.resize( M, nK);
				c.resize( nK);
				ldt.resize( nK);
				p.resize( N, nK);
				for( int k = 0, ck = 0 ; k < K ; k++){
					if( !isinf( ldt2(k)) && !isnan( ldt2(k)) ){
						learn(ck) = learn2(k);
						ldt(ck) = ldt2(k);
						c(ck) = c2(k);
						for( int i = 0 ; i < M ; i++){
							m(i,ck) = m2(i,k);
							is(i,ck) = is2(i,k);
						}
						ck++;
					}
				}
				K = nK;
			}
		}

	}

private:
	// Evaluate likelihoods on data
	void likelihoods( const array<T> &x, array<T> &p)
	{
		// Check sizes and allocate output
		int M = x.n, N = x.m;
		if( M != m.m)
			throw std::runtime_error( "gmm_t::likelihoods(): Incompatible sizes");
		p.resize( N);
		for( int i = 0 ; i < N ; i++)
			p(i) = log( 0.);

//#pragma omp parallel for
			for( int k = 0 ; k < K ; k++){
				T gc = log( c(k)) + 0.5*ldt(k) - 0.5*M*log(2*M_PI);
				for( int j = 0 ; j < N ; j++){
					T qt = 0;
					for( int i = 0 ; i < M ; i++)
						qt += is(i,k) * (x(j,i) - m(i,k)) * (x(j,i) - m(i,k));
					p(j) = lsum( p(j), gc - 0.5*qt);
				}
			}
	}

public:
	// Save the data
	void save( const std::string& filename)
	{
		using namespace std;
		ofstream f( filename.c_str(), ios::out | ios::binary);
		
		// number of gaussians
		f.write( (char*)&K, sizeof( int));

		// dimension
		f.write( (char*)&m.m, sizeof( int));

		// priors
		f.write( (char*)&c(0), K*sizeof( T));

		// means
		f.write( (char*)&m(0), m.m*K*sizeof( T));

		// inverse variances
		f.write( (char*)&is(0), is.m*K*sizeof( T));
	}

	// Load the data
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

		// priors
		c.resize( K);
		f.read( (char*)&c(0), K*sizeof( T));

		// means
		m.resize( M, K);
		f.read( (char*)&m(0), M*K*sizeof( T));

		// inverse variances
		is.resize( M, K);
		f.read( (char*)&is(0), M*K*sizeof( T));

		// Get the determinants
		ldt.resize( K);
		for( int k = 0 ; k < K ; k++){
			ldt(k) = 0;
			for( int i = 0 ; i < M ; i++)
				ldt(k) += log( is(i,k));
		}
	}

};

#endif
