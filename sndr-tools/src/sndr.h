// sndr.h -- Classes to perform sound classification tasks
// language: C++
// author  : Paris Smaragdis

#ifndef __SNDR_H
#define __SNDR_H

#include "intransp.h"
#include "aufeat.h"
#include "gmm.h"
#include "hmm.h"
#include "optparse.h"
#include "wavfile.h"

#include <list>
#include <vector>
#include <sys/time.h>

template <class T>
class AudioFeatureExtractor_t{
public:
	int b;   // Number of Mel coeffs
	int fb;  // Number of Mel filterbanks
	int hp;  // Hop size factor (1, 2, 4, etc) 
	int av;  // Feature averaging
	int sz;  // Analysis window in samples
	T tsz;   // Size of analysis window in seconds
	T flo;   // Lowest analysis frequency
	T fhi;   // Highest analysis frequency
	T thr;   // Silence threshold
	T srate; // Input sampling rate
	std::string fopts; // Audio feature options
	aufeat_t<T> F; // Feature extractor

	// Initializers
	AudioFeatureExtractor_t( int _b = 13, T _tsz = .1, T _flo = 120, T _fhi = 6855,
									int _fb = 30, std::string _fopts = std::string( "cdm"), 
									T _thr = 0, int _av = 1) : 
		b(_b), fb(_fb), hp( 1), av( _av), sz( 0), tsz(_tsz),
		flo(_flo), fhi(_fhi), thr(_thr), srate(0), fopts(_fopts) {}

	AudioFeatureExtractor_t( AudioFeatureExtractor_t &f) : 
		b(f.b), fb(f.fb), hp(f.hp), av(f.av), sz(0), tsz(f.tsz),
		flo(f.flo), fhi(f.fhi), thr(f.thr), srate(0), fopts(f.fopts), F(f.F) {}

	void report()
	{
		using namespace std;
		cout << "b = " << b << endl;
		cout << "fb = " << fb << endl;
		cout << "hp = " << hp << endl;
		cout << "av = " << av << endl;
		cout << "sz = " << sz << endl;
		cout << "tsz = " << tsz << endl;
		cout << "flo = " << flo << endl;
		cout << "fhi = " << fhi << endl;
		cout << "thr = " << thr << endl;
		cout << "srate = " << srate << endl;
		cout << "fopts = " << fopts << endl;
		F.report();
	}

	// Check options similarity
	bool operator==( const AudioFeatureExtractor_t<T> &A) const
	{
		return (b == A.b) && (fb == A.fb) && (hp == A.hp) && (av == A.av) &&
			(tsz == A.tsz) && (flo == A.flo) && (fhi == A.fhi) && (thr == A.thr) &&
			(srate == A.srate) && (fopts.compare( A.fopts) == 0);
	}
	bool operator!=( const AudioFeatureExtractor_t<T> &A) const { return !(*this == A); }
									
	// Extract sound features
	void operator()( const array<T> &s, int sr, array<T> &f, array<int> &p, bool thrm = false)
	{
	using namespace std;
		// Init the feature structure according to the new sample rate
		if( sr != srate){
			sz = pow( 2., int( log2( tsz*sr)));
			F.setup( sr, sz, round( 3*log(sr)), b, flo, fhi, fopts);
//		F.report();
			srate = sr;
		}

		// Some preallocation
		array<T> e( (s.size()-sz)/(sz/hp)+1);
		
#if 0
		// Get features of input sound using online estimation (has 4-frame delay)
		f.resize( F.o, (s.size()-sz)/(sz/hp)+1);
		for( int i = 0, j = 0 ; i < s.size()-sz ; i+=sz/hp, j++){
			array<T> st( s.v+i, sz), ft( f.v+j*F.o, F.o);
			e(j) = F.extract( st, ft);
		}
#else
		// Get features of input sound using offline estimation
		F.extract_offline( f, e, s, sz/hp);
#endif
	
		// Mark low-energy frames if threshold is non-zero
		if( thr){
			// Get peak volume
			T pk = 0;
			for( size_t i = 0 ; i < e.size() ; ++i)
				pk = max( pk, e(i));
			
			// Find passable frames
			p.resize( e.size());
			int pc = 0;
			for( size_t i = 0 ; i < e.size() ; ++i){
				p(i) = e(i) >= thr*pk;
				pc += p(i);
			}

			// Remove silent frames if proper flag is set
			if( thrm){
				// Back up the data
				array<T> f2( f.m, f.n);
				for( size_t i = 0 ; i < f.size() ; ++i)
					f2(i) = f(i);
				
				// Keep only the loud parts
				f.resize( F.o, pc);
				for( size_t i = 0, j = 0 ; i < e.size() ; ++i){
					if( p(i)){
						for( int k = 0 ; k < F.o ; ++k)
							f(k,j) = f2(k,i);
						++j;
					}
				}
				std::cout << "Volume trimmed from " << f2.n << " frames to " << f.n << " frames" << std::endl;
			}
		}
		
		// Feature averaging?
		if( av > 1){
			// Back up the data
			array<T> f2( f.m, f.n);
			for( size_t i = 0 ; i < f.size() ; ++i)
				f2(i) = f(i);
			
			// Start averaging
			f.resize( F.o, f.n/av);
			for( size_t i = 0 ; i < f.m ; ++i)
				for( size_t j = 0 ; j < f.n ; ++j){
					f(i,j) = 0;
					for( int k = 0 ; k < av ; ++k)
						if( j+k < f2.n)
							f(i,j) += f2(i,j+k);
					f(i,j) /= av;
				}
		}
		
		// Transpose in place to make the cache happy during training
		intp( &f(0), f.m, f.n);
		f.k = f.m; f.m = f.n; f.n = f.k; f.k = 1;

#if 0
		// Dump to debugging file
		{
			ofstream df( "/tmp/sndr-debugging-dump.dat", ios::out);
			df.write( ( char*)&f.m, sizeof( int));
			df.write( ( char*)&f.n, sizeof( int));
			df.write( ( char*)f.v, sizeof( T)*f.m*f.n);
		}
#endif
	}
};

// Audio feature container
template <class T>
class AudioFeatures_t{
public:
	AudioFeatureExtractor_t<T> &F; // Feature extractor to use
	array<T> C; // Consolidated features
	std::list<array<T> > D; // Feature data
	std::list<array<int> > S; // Silent frame data

	// Initialize from a feature extractor
	AudioFeatures_t( AudioFeatureExtractor_t<T> &_F) : F(_F) {}

	// Append sound to learning dictionary
	void operator()( const array<T> &in, T sr, bool thrm = false)
	{
		D.push_back( array<T>());
		S.push_back( array<int>());
		F( in, sr, D.back(), S.back(), thrm);
	}

	// Consolidate all feature sets
	void consolidate()
	{
		// Figure out how much data we have
		int fs = 0;
		for( typename std::list<array<T> >::iterator i = D.begin() ; i != D.end() ; ++i)
			fs += (*i).m;

		// Nothing to do
//		if( fs == 0)
//			throw std::runtime_error( "AudioFeatures_t<T>::consolidate(): Can't consolidate data, list is empty");

		// Consolidate features
		int ck = 0;
		const int del = F.fopts.find( 'd') != std::string::npos;
		C.resize( fs-del*5*D.size(), D.front().n);
		while( D.size()){
			for( unsigned int k = del*5 ; k < D.front().m ; ++k, ++ck)
				for( unsigned int j = 0 ; j < D.front().n ; ++j)
					C(ck,j) = D.front()(k,j);
			D.pop_front();
			S.pop_front();
		}
		std::cout << "Overall feature size " << C.m << " x " << C.n << std::endl;
	}

	// Clear buffers
	void clear()
	{
		C.resize( 0);
		while( !D.empty())
			D.pop_front();
		while( !S.empty())
			S.pop_front();
	}
};


// Single-class model
template <class T>
class AudioModel_t{
public:

#ifdef __HMM_TRAIN
	hmm_t<T> G; // HMM model for the sound class
#else
	gmm_t<T> G; // GMM model for the sound class
#endif
	AudioFeatureExtractor_t<T> &F; // Feature extractor reference, used as a test

	// Constructors
	AudioModel_t( AudioFeatureExtractor_t<T> &_F) : F(_F) {}
	AudioModel_t( AudioFeatureExtractor_t<T> &_F, const int k) : G(k), F(_F) {}

	// Train G from data A
	void operator()( AudioFeatures_t<T> &A, const int it)
	{
		if( A.F != F)
			throw std::runtime_error( "AudioModel_t trained feature space doesn't match.");
		
		// Consolidate features
		A.consolidate();

		// Learn the model
		G.train( A.C, it);
		std::cout << "Trained model in " << it << " iterations." << std::endl;
	}
	
	// Train G from data A and an initial model
#ifdef __HMM_TRAIN
	void operator()( AudioFeatures_t<T> &A, const int it, hmm_t<T> Gb)
#else
	void operator()( AudioFeatures_t<T> &A, const int it, gmm_t<T> Gb)
#endif
	{
		if( A.F != F)
			throw std::runtime_error( "AudioModel_t model-trained feature space doesn't match.");
		
		// Consolidate features
		A.consolidate();

		// Learn the model
		G.train( A.C, it, Gb);
		std::cout << "Updated model in " << it << " iterations." << std::endl;
	}
	
	// Save the model to disk
	bool save( const std::string &f)
	{
		if (f.empty()) {
		  std::cout << "error: AudioModel_t::save(\"\");" << std::endl;
		  return false;
		}

		G.save( f);

		// Sneakily append the feature data to the model file
		std::ofstream ff( f.c_str(), std::ios::binary | std::ios::app);
		ff.seekp( 0, std::ios::end);
		ff.write( "feat", 4*sizeof( char));
		ff.write( (char*)&F.b, sizeof( int));
		ff.write( (char*)&F.fb, sizeof( int));
		ff.write( (char*)&F.hp, sizeof( int));
		ff.write( (char*)&F.av, sizeof( int));
		ff.write( (char*)&F.tsz, sizeof( T));
		ff.write( (char*)&F.flo, sizeof( T));
		ff.write( (char*)&F.fhi, sizeof( T));
		ff.write( (char*)&F.thr, sizeof( T));
		ff.write( (char*)&F.srate, sizeof( T));
		ff.write( F.fopts.c_str(), F.fopts.size()*sizeof( char));

		std::cout << "Saved AudioModel_t model " << f << "." << std::endl;
		return true;
	}

	// Load a model from disk
	bool load( const std::string &f)
	{
		if (f.empty()) {
		  std::cout << "error: AudioModel_t::load(\"\");" << std::endl;
		  return false;
		}
		// Load the model
		G.load( f);
		std::cout << "Loaded model " << f << "." << std::endl;

		// Load the extra feature parameter data (if any)
		char opt[32];
		char *opp = opt;
		std::ifstream ff( f.c_str(), std::ios::binary);
		do{
			ff.read( opt, 4*sizeof( char));
			if (!ff) {
			  std::cout << "error: read only " << ff.gcount() << " of 4 chars from file '" << f << "'." << std::endl;
			  break;
			}
			ff.seekg( -3, std::ios::cur);
		}while( !ff.eof() && !(opt[0] == 'f' && opt[1] == 'e' && opt[2] == 'a' && opt[3] == 't'));
		opt[4] = '\0';
		if( strcmp( opt, "feat") == 0){
			std::cout << "Found feature parameters in AudioModel_t file " << f << "." << std::endl;
			ff.seekg( 3, std::ios::cur);
			ff.read( (char*)&F.b, sizeof( int));
			ff.read( (char*)&F.fb, sizeof( int));
			ff.read( (char*)&F.hp, sizeof( int));
			ff.read( (char*)&F.av, sizeof( int));
			ff.read( (char*)&F.tsz, sizeof( T));
			ff.read( (char*)&F.flo, sizeof( T));
			ff.read( (char*)&F.fhi, sizeof( T));
			ff.read( (char*)&F.thr, sizeof( T));
			ff.read( (char*)&F.srate, sizeof( T));
			while( !ff.eof())
				ff.read( opp++, sizeof( char));
			*--opp = '\0';
			F.fopts = opt;
		}

		// Ensure feature class initialization
		F.srate = 0;
		return true;
	}
};

// Multiple-class combiner
template <class T>
class AudioClassifier_t{
public:
	int m;  // Length of state output median filter
	T mw;   // Median filter bias
	
	T trans;       // HMM transition likelihood	
	array<T> r;    // Transition bias between states
	array<T> bias; // Likelihood bias for each HMM state

	array<int> o; // Classification output
	hmm_t<T> H;   // The master HMM
	AudioFeatureExtractor_t<T> &F; // Feature extractor reference

	// Initialize
	AudioClassifier_t( AudioFeatureExtractor_t<T> &_F) 
		: m(0), mw(1.), trans(-60), F(_F) {}

	// Combine multiple saved sound models into one classifier
	void combine( const array<std::string > &f)
	{
		// Load models in a temporary array
		std::list<AudioModel_t<T> > Al;
		for( size_t i = 0 ; i < f.size() ; ++i){
			if (f(i).empty()) {
			  std::cout << "AudioClassifier_t::combine() skipping empty filename." << std::endl;
			  continue;
			}
			Al.push_back( AudioModel_t<T>( F));
			if (!Al.back().load( f(i))) {
			  std::cout << "AudioClassifier_t::combine() failed to load model '" << f(i) << "'." << std::endl;
			}
		}

		// Combine them
		combine( Al);
	}

	// Combine multiple sound models into one classifier
	void combine( std::list<AudioModel_t<T> > &Al)
	{
#ifdef __HMM_TRAIN
		// Pack into a single HMM.
		if (Al.size() != 2)
		  throw std::runtime_error("BUG: hmm_t::combine() ignores intermediate elements of list of HMMs.");
		// todo: use all elements of list.
		::combine( H, Al.front().G, Al.back().G, 1.-trans, 1.-trans);
#else
		// Pack all GMMs in a HMM
		int M = -1;
		int ai = 0;
		for( typename std::list<AudioModel_t<T> >::iterator A = Al.begin() ; A != Al.end() ; ++A, ai++){
			// Make sure all models are relevant
//			if( (*A).F != F)
//				throw std::runtime_error( "AudioClassifier_t<T>::combine(): Input models are not using the same features");

			// Get the GMM inside the model
			gmm_t<T> &G = (*A).G;

			// Make sure that we can fit all that stuff
			if( ai == 0){
				M = G.m.m;
				H.S = Al.size();
				H.K = G.K;
				std::cout << "Making a " << H.S << "-state HMM with " << H.K << " gaussians per state." << std::endl;
				H.lPi.resize( H.S);
				H.lA.resize( H.S, H.S);
				H.c.resize( H.K, H.S);
				H.ldt.resize( H.K, H.S);
				H.m.resize( M, H.K, H.S);
				H.is.resize( M, H.K, H.S);
			}

			// Copy over the GMM to current state
			for( int k = 0 ; k < G.K ; k++){
				H.c(k,ai) = G.c(k);
				H.ldt(k,ai) = G.ldt(k);
				for( int i = 0 ; i < M ; i++){
					H.m(i,k,ai) = G.m(i,k);
					H.is(i,k,ai) = G.is(i,k);
				}
			}

			// Make the transition matrix row
			if( trans > 0){
				for( int i = 0 ; i < H.S ; i++)
					H.lA(ai,i) = log( (1.-trans)/(H.S-1));
				H.lA(ai,ai) = log( trans);
			}else{
				for( int i = 0 ; i < H.S ; i++)
					H.lA(ai,i) = trans;
				H.lA(ai,ai) = log( .9999); // **** improper value, should be ok further down when I normalize lA
			}			
		}

		// Make initial probabilities
		for( int i = 0 ; i < H.S ; i++)
			H.lPi(i) = log( 1./H.S);
		
		// Norm the priors
		for( int i = 0 ; i < H.S ; i++){
			T s = 0;
			for( int k = 0 ; k < H.K ; k++)
				s += H.c(k,i);
			for( int k = 0 ; k < H.K ; k++)
				H.c(k,i) -= log( s);
		}

		// Bias the transitions
		if( int(r.size()) == H.S*H.S){
			std::cout << "Biasing" << std::endl;
			for( int i = 0 ; i < H.S*H.S ; i++)
				H.lA(i) -= r(i);
		}

		// Normalize the transition matrix again
		for( int i = 0 ; i < H.S ; i++){
			T ls = log( 0.);
			for( int j = 0 ; j < H.S ; j++)
				ls = std::max( ls, H.lA(i,j)) + log1p( exp( -fabs( ls-H.lA(i,j))));
			for( int j = 0 ; j < H.S ; j++)
				H.lA(i,j) -= ls;
		}
#endif
		std::cout << "Transition matrix is:\n" << H.lA;
	}

	// Classify a new input
	void operator()( const array<T> &in, int sr)
	{
		// Get the sound features
		array<T> D;
		array<int> S;
//		F.report();
		F( in, sr, D, S, false);

		// Classify
		H.classify( D, o, bias);
		std::cout << "Input is " << in.size() << " window is " << F.sz << " out is " << o.size() << std::endl;

		// Relabel the silent parts to the last known class
		if( F.thr > 0){
			int c = 0;
			for( size_t i = 0 ; i < o.size() ; ++i){
				if( S(i))
					c = o(i);
				else
					o(i) = c;
			}
		}

		// Pass through a median filter to smooth out sustained sections
		if( m > 1){
			array<int> o2( o.size());
			for( size_t i = 0 ; i < o.size() ; ++i)
				o2(i) = o(i);
			int c[2];
			for( size_t i = m ; i < o.size()-m ; ++i){
				c[0] = c[1] = 0;
				for( int j = -m ; j <= m ; ++j)
					if( o2(i+j) != -1)
						++c[o2(i+j)];
				o(i) = mw*c[0] > c[1] ? 0 : 1;
			}
		}

		// Show the rates
		array<int> hs( H.S);
		for( int j = 0 ; j != H.S ; ++j)
			hs(j) = 0;
		for( size_t j = 0 ; j != o.size() ; ++j)
			hs(o(j))++;
		std::cout << "Class result histogram: " << hs << std::endl;
	}

	// Set state filtering parameters
	void filter( int _m, T _w)
	{
		m = _m;
		mw = _w;
	}

	// Make an EDL file for class ii
	void make_edl( const std::string &edl, int ii = 0)
	{
		std::ofstream e( edl.c_str());
		bool cl = false;
		if( o(0) == ii){
			e << 0 << ' ';
			cl = true;
		}
		for( size_t i = 1 ; i < o.size() ; ++i)
			if( o(i) != o(i-1)){
				e << T(i*F.sz/F.hp)/F.srate;
				if( !cl && o(i) == ii){
					e << ' ';
					cl = true;
				}else{
					e << " 0\n";
					cl = false;
				}
			}
		if( cl)
			e << T(o.size()*F.sz/F.hp)/F.srate << " 0\n";
	}

	// Make a soundfile for class ii
	void make_snd( const array<T> &x, const std::string &f, int ii = 0)
	{
		wavfile_t sf( f, std::ios::out, 1, F.srate);
		int N = 0;
		for( size_t k = 0 ; k < o.size() ; ++k)
			if( o(k) == ii){
				array<T> tt( x.v+k*F.sz/F.hp, F.sz/F.hp);
				sf.write( tt);
				++N;
			}
		std::cout << "Dumped " << N << " of " << o.size() << " frames to " << f << "." << std::endl;
	}

	// Load preexisting HMM model (potentially with feature info)
	void load( const std::string &modin)
	{
		// Load the HMM (throw exception on failure).
		H.load( modin);

		// Load any extra feature parameters.
		char opt[32];
		char *opp = opt;
		std::ifstream ff( modin.c_str(), std::ios::binary);
		do{
			ff.read( opt, 4*sizeof( char));
			ff.seekg( -3, std::ios::cur);
		}while( !ff.eof() && !(opt[0] == 'f' && opt[1] == 'e' && opt[2] == 'a' && opt[3] == 't'));
		opt[4] = '\0';
		if( strcmp( opt, "feat") == 0){
			std::cout << "HMM file " << modin << " has feature parameters." << std::endl;
			ff.seekg( 3, std::ios::cur);
			ff.read( (char*)&F.b, sizeof( int));
			ff.read( (char*)&F.fb, sizeof( int));
			ff.read( (char*)&F.hp, sizeof( int));
			ff.read( (char*)&F.av, sizeof( int));
			ff.read( (char*)&F.tsz, sizeof( T));
			ff.read( (char*)&F.flo, sizeof( T));
			ff.read( (char*)&F.fhi, sizeof( T));
			ff.read( (char*)&F.thr, sizeof( T));
			ff.read( (char*)&F.srate, sizeof( T));
#if 0
			int tbs;
			ff.read( (char*)&tbs, sizeof( int));
			if( tbs > 0){
				bias.resize( tbs);
				ff.read( (char*)bias.v, bias.size()*sizeof( T));
			}
#endif
			while( !ff.eof())
				ff.read( opp++, sizeof( char));
			*--opp = '\0';
			F.fopts = opt;
		}

		// Ensure feature class initialization.
		F.srate = 0;
	}

	// Combine multiple sound model files into an HMM model file (and pack in the feature info as well)
	void combine_models( const array<std::string> &modin, const std::string &modout)
	{
		// Load GMM models and combine them into one HMM
		combine( modin);

		// Save the HMM (throw exception on failure).
		H.save( modout);

		// Be sneaky and append the feature data on the hmm file
		std::ofstream ff( modout.c_str(), std::ios::binary | std::ios::app);
		ff.seekp( 0, std::ios::end);
		ff.write( "feat", 4*sizeof( char));
		ff.write( (char*)&F.b, sizeof( int));
		ff.write( (char*)&F.fb, sizeof( int));
		ff.write( (char*)&F.hp, sizeof( int));
		ff.write( (char*)&F.av, sizeof( int));
		ff.write( (char*)&F.tsz, sizeof( T));
		ff.write( (char*)&F.flo, sizeof( T));
		ff.write( (char*)&F.fhi, sizeof( T));
		ff.write( (char*)&F.thr, sizeof( T));
		ff.write( (char*)&F.srate, sizeof( T));
//		int bs = F.bias.size();
//		ff.write( (char*)&F.bs, sizeof( int));
//		ff.write( (char*)F.bias.v, F.bias.size()*sizeof( T));
		ff.write( F.fopts.c_str(), F.fopts.size()*sizeof( char));
	}
};

#endif
