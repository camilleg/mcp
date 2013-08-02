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
#include "aquat.h"


//
// Hidden Markov Model class
//

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

	// log sum
	inline T lsum( T x, T y)
	{
		using namespace std;
		if( fabs( x-y) > 30)
			return max( x, y);
		if( x == y)
			return x+log(2.);
		return max( x, y) + log1p( exp( -fabs( x-y)));
	}

	// Learn data
	void train( const array<T> &x, int iters = 100, const hmm_t<T> &H = hmm_t<T>())
	{
		using namespace std;

		// Input dims are reversed
		int M = x.n, N = x.m;
		cout << "input is " << M << 'x' << N << endl;

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
			for( int s = 0 ; s < S ; s++){
				for( int k = 0 ; k < K ; k++){
					ldt(k,s) = 0;
					c(k,s) = 1./K;
					for( int i = 0 ; i < M ; i++){
						m(i,k,s) = T( rand())/RAND_MAX - .5;
						is(i,k,s) = .1;
						ldt(k,s) += log( is(i,k,s));
					}
				}
			}

			// Initial values for initial and transition probabilities
			for( int s = 0 ; s < S ; s++)
				lPi(s) = log( 1./S);
			for( int i = 0 ; i < S ; i++)
				for( int j = 0 ; j < S ; j++)
					lA(i,j) = log( T(rand())/RAND_MAX);
			for( int i = 0 ; i < S ; i++){
				T ls = lA(i,0);
				for( int j = 0 ; j < S ; j++)
					ls = lsum( ls, lA(i,j));
				for( int j = 0 ; j < S ; j++)
						lA(i,j) = lA(i,j) - ls;
			}
		}else{
			
			// Copy values for Gaussians
			for( int s = 0 ; s < S ; s++){
				for( int k = 0 ; k < K ; k++){
					ldt(k,s) = H.ldt(k,s);
					c(k,s) = H.c(k,s);
					for( int i = 0 ; i < M ; i++){
						m(i,k,s) = H.m(i,k,s);
						is(i,k,s) = H.is(i,k,s);
					}
				}
			}

			// Copy values for initial and transition probabilities
			for( int i = 0 ; i < S ; i++){
				lPi(i) = H.lPi(i);
				for( int j = 0 ; j < S ; j++)
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

		// Start iterating
		for( int it = 0 ; it < iters ; it++){
		
			//
			// Expectation step
			//

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
					T tp = log( 0.);
					for( int k = 0 ; k < K ; k++)
						tp = lsum( tp, q(j,k,s));
					lp(s,j) = tp;
				}

			// Get alphas
			for( int i = 0 ; i < S ; i++)
				la(i,0) = lPi(i) + lp(i,0);
			for( int t = 0 ; t < N-1 ; t++)
				for( int j = 0 ; j < S ; j++){
					T ls = log( 0.);
					for( int i = 0 ; i < S ; i++)
						ls = lsum( ls, la(i,t) + lA(i,j));
					la(j,t+1) = ls + lp(j,t+1);
				}

			// Get betas
			for( int i = 0 ; i < S ; i++)
				lb(i,N-1) = log( 0.);
			lb(S-1,N-1) = 0;
			for( int t = N-2 ; t > -1 ; t--)
				for( int i = 0 ; i < S ; i++){
					T ls = log( 0.);
					for( int j = 0 ; j < S ; j++)
						ls = lsum( ls, lb(j,t+1) + lA(i,j) + lp(j,t+1));
					lb(i,t) = ls;
				}

			// Get Xi
			for( int i = 0 ; i < S*S ; i++)
				xi(i) = log( 0.);
			for( int t = 0 ; t < N-1 ; t++){
				T ls = log( 0.);
				for( int i = 0 ; i < S ; i++)
					for( int j = 0 ; j < S ; j++){
						txi(i,j) = lA(i,j) + la(i,t) + lp(j,t+1) + lb(j,t+1);
						ls = lsum( ls, txi(i,j));
					}
				for( int i = 0 ; i < S*S ; i++)
					xi(i) = lsum( xi(i), txi(i) - ls);
			}

			// Get gamma
			for( int j = 0 ; j < N ; j++){
				T ls = log( 0.);
				for( int s = 0 ; s < S ; s++)
					ls = lsum( ls, la(s,j)+lb(s,j));
				for( int s = 0 ; s < S ; s++){
					T tg = la(s,j) + lb(s,j) - lp(s,j) - ls;
					for( int k = 0 ; k < K ; k++)
						g(j,k,s) = tg + q(j,k,s);
				}
			}

			// Get overall likelihood
			lk(it) = log( 0.);
			for( int i = 0 ; i < S ; i++)
				lk(it) = lsum( lk(it), la(i,N-1));
			if( 1){ //( !((it+1)%25) | it == iters-1){
				cout << "Iteration: " << it+1 << " Likelihood: " << lk(it) << endl;
				if( it == 0)
					aq_window( rand(), "HMM Training Likelihood");
				aq_plot( &lk(0), it);
				char ts[64]; sprintf( ts, "Iteration: %d, Likelihood: %.2f", it, lk(it));
				aq_text( ts, .5, 1./15);
			}

			// Get out of log domain
			for( int i = 0 ; i < N*K*S ; i++)
				g(i) = exp( g(i));


			//
			// Maximization step
			//

			// Initial probabilities
			for( int s = 0 ; s < S ; s++){
				T tp = log( .0);
				for( int i = 0 ; i < K ; i++)
					tp = lsum( tp, log( g(0,i,s)));
				lPi(s) = tp;
			}

			// Transition matrix
			for( int i = 0 ; i < S ; i++){
				T ls = log( 0.);
				for( int j = 0 ; j < S ; j++)
					ls = lsum( ls, xi(i,j));
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

	// Classify based on an known HMM
	void classify( const array<T> &x, array<int> &q, array<T> &bias = array<T>(), int ist = -1)
	{
		using namespace std;
		int M = x.n, N = x.m;

		// Get state probabilities
		array<T> lB( S, N);
		for( int i = 0 ; i < S*N ; i++)
			lB(i) = -HUGE_VAL;

#pragma omp parallel for
		for( int s = 0 ; s < S ; s++)
			for( int k = 0 ; k < K ; k++){
				T gc = log( c(k,s)) + 0.5*ldt(k,s) - 0.5*M*log(2*M_PI);
				for( int j = 0 ; j < N ; j++){
					T qt = 0;
					for( int i = 0 ; i < M ; i++)
						qt += is(i,k,s) * (x(j,i) - m(i,k,s)) * (x(j,i) - m(i,k,s));
					lB(s,j) = lsum( lB(s,j), gc - 0.5*qt);
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

		// Show likelihoods
		//aq_window( rand(), "HMM State Posteriors"); aq_mplot( &lB(0), lB.m, lB.n);

		// Run viterbi
		viterbi( lB, q, ist);
	}

	// Viterbi decoding
	void viterbi( const array<T> &lB, array<int> &q, int ist = -1)
	{
		using namespace std;
		int N = lB.n;

		// Temp space
		array<T> d( S, 2);
		array<int> p( S, N);

		// Override initial state if one is provided
		array<double> nlPi(S);
		if( ist != -1){
			for( int i = 0 ; i < S ; i++)
				nlPi(i) = log( 0.);
			nlPi(ist) = 0.;
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
		using namespace std;
		int M = lB.m, N = lB.n;

		// Allocate output
		q.resize( N);

		// Initial probability vector
		array<T> pp( S);
		for( int i = 0 ; i < S ; i++)
			pp(i) = lPi(i);

		// Start going at it
		int a = 0, b = 1;
		while( b < N){

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

			// If there was an advance keep it, otherwise increase window and repeat
			if( t > 0){
				for( int i = 0 ; i < t ; i++)
					q(a+i) = s(0,i);
				for( int i = 0 ; i < S ; i++)
					pp(i) = lA(q(a+t-1),i);
				a += t;
				b = a+1;
			}else
				b++;
		}
	}

	// Interim Viterbi decoding
	void iviterbi( const array<T> &lB, int a, int b, array<T> &lP, array<int> &q)
	{
		using namespace std;
		int N = b-a;

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

	// Save the model
	void save( const std::string& filename)
	{
		using namespace std;
		ofstream f( filename.c_str(), ios::out | ios::binary);
		
		// number of states
		f.write( (char*)&S, sizeof( int));

		// initial log probabilities
		f.write( (char*)&lPi(0), S*sizeof( T));

		// log transition matrix
		f.write( (char*)&lA(0), S*S*sizeof( T));

		// the number of gaussians
		f.write( (char*)&K, sizeof( int));

		// dimension
		int M = m.size()/(K*S);
		f.write( (char*)&M, sizeof( int));

		// priors
		f.write( (char*)&c(0), K*S*sizeof( T));

		// means
		f.write( (char*)&m(0), M*K*S*sizeof( T));

		// inverse variances
		f.write( (char*)&is(0), M*K*S*sizeof( T));

		cout << "Saved hmm_t model '" << filename << "'.\n";
	}

	// Load a model
	bool load( const std::string& filename)
	{
		using namespace std;
		ifstream f( filename.c_str(), ios::in | ios::binary);
		// todo: error handling
		
		// number of states
		f.read( (char*)&S, sizeof( int));
		if (S <= 0) {
		  cout << "Problem loading HMM File '" << filename << "'.\n";
		  throw std::runtime_error( "load(): nonpositive number of states");
		  return false;
		}

		// initial log probabilities
		lPi.resize( S);
		f.read( (char*)&lPi(0), S*sizeof( T));

		// log transition matrix
		lA.resize( S, S);
		f.read( (char*)&lA(0), S*S*sizeof( T));

		// the number of gaussians
		f.read( (char*)&K, sizeof( int));
		if (K <= 0) {
		  cout << "Problem loading HMM File '" << filename << "'.\n";
		  throw std::runtime_error( "load(): nonpositive number of gaussians");
		  return false;
		}

		// dimension
		const int M = m.size()/(K*S);
		if (M <= 0) {
		  cout << "Problem loading HMM File '" << filename << "'.\n";
		  throw std::runtime_error( "load(): nonpositive number of dimensions");
		  return false;
		}
		f.read( (char*)&M, sizeof( int));

		// priors
		c.resize( K, S);
		f.read( (char*)&c(0), K*S*sizeof( T));

		// means
		m.resize( M, K, S);
		f.read( (char*)&m(0), M*K*S*sizeof( T));

		// inverse variances
		is.resize( M, K, S);
		f.read( (char*)&is(0), M*K*S*sizeof( T));

		// Compute the determinants.
		ldt.resize( K, S);
		for( int s = 0 ; s < S ; s++)
			for( int k = 0 ; k < K ; k++){
				ldt(k,s) = 0;
				for( int i = 0 ; i < M ; i++)
					ldt(k,s) += log( is(i,k,s));
			}
		return true;
	}

};

// Combine two HMMs in a loop
template <class T>
void combine( hmm_t<T> &H, const hmm_t<T> &h1, const hmm_t<T> &h2, T p1 = .5, T p2 = .5)
{
	using namespace std;
	H.S = h1.S + h2.S;
	H.K = h1.K;
	if( H.K != h2.K)
		throw std::runtime_error( "combine(): HMM state K's are incompatible");

	// Allocate the model parameters
	int M = h1.m.m;
	if( M != h2.m.m)
		throw std::runtime_error( "combine(): Input sizes are incompatible");
	H.lPi.resize( H.S);
	H.lA.resize( H.S, H.S);
	H.c.resize( H.K, H.S);
	H.ldt.resize( H.K, H.S);
	H.m.resize( M, H.K, H.S);
	H.is.resize( M, H.K, H.S);

	// Copy over the model data
	for( int s = 0 ; s < h1.S ; s++)
		for( int k = 0 ; k < h1.K ; k++){
			H.c(k,s) = h1.c(k,s);
			H.ldt(k,s) = h1.ldt(k,s);
			for( int i = 0 ; i < M ; i++){
				H.m(i,k,s) = h1.m(i,k,s);
				H.is(i,k,s) = h1.is(i,k,s);
			}
		}
	for( int s = 0 ; s < h2.S ; s++)
		for( int k = 0 ; k < h2.K ; k++){
			H.c(k,h1.S+s) = h2.c(k,s);
			H.ldt(k,h1.S+s) = h2.ldt(k,s);
			for( int i = 0 ; i < M ; i++){
				H.m(i,k,h1.S+s) = h2.m(i,k,s);
				H.is(i,k,h1.S+s) = h2.is(i,k,s);
			}
		}
	
	// Make transition matrix and initial probabilities
	for( int i = 0 ; i < H.lA.size() ; i++)
		H.lA(i) = -HUGE_VAL;
	for( int i = 0 ; i < h1.S ; i++){
		H.lPi(i) = h1.lPi(i);
		for( int j = 0 ; j < h1.S ; j++)
			H.lA(i,j) = h1.lA(i,j);
	}
	for( int i = 0 ; i < h2.S ; i++){
		H.lPi(h1.S+i) = h2.lPi(i);
		for( int j = 0 ; j < h2.S ; j++)
			H.lA(h1.S+i,h1.S+j) = h2.lA(i,j);
	}
	H.lA(h1.S-1,h1.S) = log( p1);
	H.lA(h1.S+h2.S-1,0) = log( p2);

	// Normalize them
	T ps = 0;
	for( int i = 0 ; i < H.S ; i++)
		ps += exp( H.lPi(i));
	for( int i = 0 ; i < H.S ; i++)
		H.lPi(i) = log( exp( H.lPi(i))/ps);

	for( int i = 0 ; i < h1.S ; i++)
		H.lA(h1.S-1,i) += log( (1.-p1));
	for( int i = 0 ; i < h2.S ; i++)
		H.lA(h1.S+h2.S-1,h1.S+i) += log( (1.-p2));
}

#endif
