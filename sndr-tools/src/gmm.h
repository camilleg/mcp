// gmm.h -- Gaussian Mixture Model class
// language: C++
// author  : Paris Smaragdis

#ifndef __GMM_H__
#define __GMM_H__

#include <cmath>
#include <fstream>
#include <algorithm>
#include <iostream>

#include "pmmintrin.h"

#include "array.h"
#include "aquat.h"

// Gaussian mixture model class
template <class T>
class gmm_t {
public:

	int K; // Gaussians
	T dg; // Diagonal load
	array<T> ldt, c, m, is;

	// Default values
	gmm_t( int k = 0, T d = __FLT_EPSILON__) : K( k), dg( d) {}

	// log sum
	inline T lsum( T x, T y)
	{
		using namespace std;
		if( x == y)
			return x+log(2.);
		return max( x, y) + log1p( exp( -fabs( x-y)));
	}

	// Learn data
	void train( const array<T> &x, int iters = 100, 
						  const gmm_t<T> &G = gmm_t<T>(), bool prior = false)
	{
		using namespace std;

		// Remember the sizes
		int M = x.n, N = x.m;

		// Do a check on the input
		for( int i = 0 ; i < M*N ; i++)
			if( isinf( x(i)) | isnan( x(i)))
				throw std::runtime_error( "gmm_t::train(): Found NaN/Inf in input");

		// Setup arrays
		ldt.resize( K);
		c.resize( K);
		m.resize( M, K);
		is.resize( M, K);
		array<T> p( N, K);
		array<T> lk( iters);

		// Randomize
#ifndef __NO_RAND_SEED
		time_t tm; time( &tm);
		srand( tm);
		(void)rand();
#endif

		// Sort out the learning situation
		array<int> learn( K);
		if( G.K & prior){
			for( int k = 0 ; k < K ; k++)
				learn(k) = k >= G.K;
		}else
			for( int k = 0 ; k < K ; k++)
				learn(k) = 1;

		// Initial values
		for( int k = 0 ; k < K ; k++){
			if( learn( k)){
				ldt(k) = 0;
				c(k) = 1.;
				for( int i = 0 ; i < M ; i++){
//					m(i,k) = x(k,i);
//					m(i,k) = T( rand())/RAND_MAX - .5;
					int ri = (N-1)*double( rand())/RAND_MAX;
//					std::cout << ri << std::endl;
					m(i,k) = x(ri,i);
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
					mx = max( mx, p(j,i));
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
				cout << "Iteration: " << it+1 << " Likelihood: " << lk(it) << endl;
				if( it == 0)
					aq_window( rand(), "GMM Training Likelihood");
				aq_plot( &lk(1), it-1);
				char ts[64]; sprintf( ts, "Iteration: %d, Likelihood: %.2f", it+1, lk(it));
				aq_text( ts, .5, 1./15);
				if( it == iters-1)
					aq_close();
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
#if 1
						T ms = 0;
						for( int j = 0 ; j < N ; j++)
							ms += x(j,i) * p(j,k);
						m(i,k) = ms/ps;
#else
						m(i,k) = cblas_ddot( N, &x.v[x.m*i], 1, &p.v[p.m*k], 1)/ps;
#endif
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

			// Remove blown up states
			int nK = K;
			for( int k = 0 ; k < K ; k++)
				if( isinf( ldt(k)) | isnan( ldt(k)) )
					nK--;
			if( nK != K){
				cout << "Iteration " << it << ", " << K-nK << " blown states" << endl;
				array<T> m2( m), is2( is), c2( c), ldt2( ldt);
				array<int> learn2( learn);
				m.resize( M, nK);
				is.resize( M, nK);
				c.resize( nK);
				ldt.resize( nK);
				p.resize( N, nK);
				for( int k = 0, ck = 0 ; k < K ; k++){
					if( !isinf( ldt2(k)) & !isnan( ldt2(k)) ){
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

	// Evaluate likelihoods on data
	void likelihoods( const array<T> &x, array<T> &p)
	{
		using namespace std;
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
/*
		// Show me
		aq_window();
		aq_plot( &p(0), N);
		char ts[64];
		T mn = 0, mi = p(0), ma = p(0);
		for( int i = 0 ; i < N ; i++){
			mi = min( mi, p(i));
			ma = max( ma, p(i));
			mn += p(i);
		}
		mn /= N;
		sprintf( ts, "Min/Max/Mean: %.2f/%.2f/%.2f", mi, ma, mn);
		cout << ts << endl;
		aq_text( ts, .5, 1./15);
		aq_close();*/
	}

	// Save the data
	void save( std::string fn)
	{
		using namespace std;
		ofstream f( fn.c_str(), ios::out | ios::binary);
		
		// Write the number of gaussians
		f.write( (char*)&K, sizeof( int));

		// Write dimension
		f.write( (char*)&m.m, sizeof( int));

		// Write priors
		f.write( (char*)&c(0), K*sizeof( T));

		// Write means
		f.write( (char*)&m(0), m.m*K*sizeof( T));

		// Write inverse variances
		f.write( (char*)&is(0), is.m*K*sizeof( T));
	}

	// Load the data
	void load( std::string fn)
	{
		using namespace std;
		ifstream f( fn.c_str(), ios::in | ios::binary);
		
		// Read the number of gaussians
		f.read( (char*)&K, sizeof( int));

		// Read dimension
		int M;
		f.read( (char*)&M, sizeof( int));

		// Read priors
		c.resize( K);
		f.read( (char*)&c(0), K*sizeof( T));

		// Read means
		m.resize( M, K);
		f.read( (char*)&m(0), M*K*sizeof( T));

		// Read inverse variances
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
