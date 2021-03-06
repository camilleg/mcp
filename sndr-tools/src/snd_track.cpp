// snd_track.cpp -- Track an example sound from an audio scene
// language: C++
// author  : Paris Smaragdis

#include "intransp.h"
#include "aufeat.h"
#include "hmm.h"
#include "optparse.h"
#include "wavfile.h"
#include "real_t.h"

#include <list>
#include <sys/time.h>

using namespace std;

// Feature extraction

int fextract( const array<real_t> &s, array<real_t> &f, array<int> &p, int sr, real_t tsz, int hp, 
	int b, string opt, real_t mnf, real_t mxf, int av, real_t thr, real_t thrm = true)
{
	// Setup feature object
	aufeat_t<real_t> F;
	int sz = pow( 2., int( log2( tsz*sr)));
	F.setup( sr, sz, 3*log(sr), b, mnf, mxf, opt);
	array<real_t> e( s.size()/sz);

#if 0
	// Get features of input sound
	// real-time algorithm has a delay of 4 frames when computing deltas
	f.resize( F.o, s.size()/sz);
	for( int i = 0, j = 0 ; i < s.size()-sz ; i+=sz, ++j){
		array<real_t> st( s.v+i, sz), ft( f.v+j*F.o, F.o);
		e(j) = F.extract( st, ft);
	}
#else
	// Offline feature computation is in-time with input
	F.extract_offline( f, e, s, sz/hp);
#endif

	// Remove low volume samples, thr is fraction of the peak
	if( thr){
		// Get peak volume
		real_t pk = 0;
		for( size_t i = 0 ; i < e.size() ; ++i)
			pk = max( pk, e(i));

		// Find passable frames
		p.resize( e.size());
		int pc = 0;
		for( size_t i = 0 ; i < e.size() ; ++i){
			p(i) = e(i) >= thr*pk;
			pc += p(i);
		}

		// Remove silent frames?
		if( thrm){
			// Back up the data
			array<real_t> f2( f.m, f.n);
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
			cout << "Volume trimmed from " << f2.n << " frames to " << f.n << " frames" << endl;
		}
	}

	// Feature averaging?
	if( av > 1){
		// Back up the data
		array<real_t> f2( f.m, f.n);
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
	intp( &f[0], f.m, f.n);
	f.k = f.m; f.m = f.n; f.n = f.k; f.k = 1;
	
	// Return window size in samples
	return sz;
}

// Learning

int learn( const array<real_t> &in, const array<real_t> &s, const int K, const int it, array<real_t> &t, hmm_t &H)
{
#ifdef __HMM_TRAIN
	// Learn overall input
	hmm_t H1( K, 1);
	H1.train( in, it);
	cout << "Learned UBM" << endl;
	
	// Learn target sound
	hmm_t H2( K, 1);
	H2.train( s, it, H1);
	cout << "Learned target model" << endl;
	
	// Pack into a single HMM
	combine( H, H1, H2, real_t(1.-t(0)), real_t(1.-t(t.size()==1 ? 0 : 1)));
	cout << "Packed models in an HMM" << endl;

	// Define in and out state threshhold
	return H1.S;
#else
	gmm_t G[2] = {K,K};
	// http://stackoverflow.com/questions/18854702/instantiate-c-templated-class-with-arg-and-as-array

	G[0].train(in, it);
	cout << "Learned UBM." << endl; // overall input

	G[1].train(s, it, G[0]);
	cout << "Learned target model." << endl; // target sound

	// Make room in HMM
	H.K = G[0].numGaussians();
	H.S = 2;
	H.lPi.resize(H.S);
	H.lA.resize(H.S, H.S);
	H.smps.resize(H.S);

	// Copy GMMs into HMM
	for (int s=0; s<H.S; ++s) H.smps(s).copyGuts(G[s]);

	// Make the transition matrix row
	if( t.size() == 0)
		t.push_back( .9);
	H.lA(0,0) = log(t(0));
	H.lA(0,1) = log1p(-t(0));
	const real_t tt = t.size() > 1 ? t(1) : t(0);
	H.lA(1,0) = log1p(-tt);
	H.lA(1,1) = log(tt);
	cout << "Transition matrix: \n" << H.lA;

	// Make initial probabilities
	H.lPi(0) = log(0.5);
	H.lPi(1) = log(0.5);

	// Normalize the priors
	for (int s=0; s<H.S; ++s) H.smps(s).normalizePriors();
	cout << "Packed models in an HMM." << endl;
	return 0;
#endif
}

// Find the target

array<int> search( const array<real_t>& in, hmm_t& H, const array<real_t>& b)
{
	// Run the classification
	array<int> o;
	H.classify( o, in, b);
	return o;
}

int main( int argc, const char **argv)
{
	struct Options {
	  int b, K, it, m, av;
	  real_t nfn, mn_freq, mx_freq, tsz;
	  real_t tgs, tge, thr, hp, mw;
	  array<real_t> tr, B;
	  string sf_in, sf_targ, opt, fo, edl;
	};

	// Get options
	Options O;
	O.sf_in = getoption<string>( "-i", argc, argv, "", "Input sound"); // Input sound
	O.sf_targ = getoption<string>( "-t", argc, argv, "", "Target sound"); // Target sound
	O.tsz = getoption<real_t>( "-s", argc, argv, .05, "Window size"); // Window size (sec)
	O.hp = getoption<real_t>( "-H", argc, argv, 1, "Hop size"); // Hop size (window size/hp)
	O.tgs = getoption<real_t>( "-S", argc, argv, -1, "Target start time"); // Target start time (sec)
	O.tge = getoption<real_t>( "-E", argc, argv, -1, "Target end time"); // Target end time (sec)
	O.tr = mgetoption<real_t>( "-p", argc, argv, "State transition probabilities"); // State transition probabilities
	O.K = getoption<int>( "-K", argc, argv, 8, "Number of Gaussians"); // Gaussians to use
	O.it = getoption<int>( "-e", argc, argv, 50, "Training iterations"); // Training iterations
	O.nfn = getoption<int>( "-n", argc, argv, 30, "Number of filterbank filters"); // // Number of filters in filterbank, n*log( fs)
	O.b = getoption<int>( "-b", argc, argv, 13, "Number of coefficients"); // Coeffs to use
	O.m = getoption<int>( "-m", argc, argv, 0, "Length of state output filter"); // Length of state output filter
	O.mw = getoption<real_t>( "-w", argc, argv, 1, "Bias of state output filter"); // Bias of state output filter
	O.mn_freq = getoption<real_t>( "-l", argc, argv, 133, "Lowest frequency"); // Lowest freq
	O.mx_freq = getoption<real_t>( "-h", argc, argv, 6855, "Highest frequency"); // Highest freq
	O.B = mgetoption<real_t>( "-B", argc, argv, "Bias"); // Bias
	O.thr = getoption<real_t>( "-T", argc, argv, 0, "Peak threshold"); // Peak threshhold
	O.opt = getoption<string>( "-F", argc, argv, "mcd", "Feature options"); // Feature options
	O.fo = getoption<string>( "-f", argc, argv, "output.wav", "Output filename"); // Output name
	O.edl = getoption<string>( "-D", argc, argv, "", "EDL filename"); // EDL file name
	O.av = getoption<int>( "-a", argc, argv, 1, "Feature averaging"); // Feature averaging

	// Load the input sound
	wavfile_t sfi( O.sf_in, ios::in);
	array<real_t> in( sfi.frames);
	sfi.read_mono( in);
	cout << argv[0] << ": loaded input soundfile of " << double( sfi.frames)/sfi.samplerate << " seconds." << endl;

	// Is the target sound a subset of the input?
	array<real_t> targ;
	if( O.tgs != -1 && O.tge != -1){
		int is = sfi.samplerate*O.tgs, ie = sfi.samplerate*O.tge;
		targ.resize( ie-is);
		for( int i = 0 ; i < ie-is ; ++i)
			targ(i) = in(is+i);
		cout << argv[0] << ": located target." << endl;
	}else{
		// Load the target sound
		wavfile_t sfi2( O.sf_targ, ios::in);
		targ.resize( sfi2.frames);
		sfi2.read_mono( targ);
		cout << argv[0] << ": loaded target soundfile of " << double( sfi2.frames)/sfi2.samplerate << " seconds." << endl;

		// Check to make sure that the two sample rates are the same
		if( sfi.samplerate != sfi2.samplerate)
			throw std::runtime_error( "Sample rates of inputs are not compatible");
	}

	// Normalize the inputs
	targ.normalize_max1();
	in.normalize_max1();

	// Start tracking performance
	struct timeval time;
	gettimeofday( &time, NULL);
	double start_time = time.tv_sec + double( time.tv_usec)/1000000.;

	// Get the sound features
	array<real_t> fi, ft;
	array<int> sl;
	int sz = fextract( targ, ft, sl, sfi.samplerate, O.tsz, O.hp, O.b, O.opt, O.mn_freq, O.mx_freq, 1, O.thr, true);
	cout << argv[0] << ": extracted target features, " << ft.m << " x " << ft.n << "." << endl;
	cout << argv[0] << ": analysis window has " << sz << " samples." << endl;
	fextract( in, fi, sl, sfi.samplerate, O.tsz, O.hp, O.b, O.opt, O.mn_freq, O.mx_freq, O.av, O.thr, false);
	cout << argv[0] << ": extracted input features, " << fi.m << " x " << fi.n << "." << endl;

	// Do the learning
	hmm_t H;
	int stth = learn( fi, ft, O.K, O.it, O.tr, H);
	cout << argv[0] << ": learned." << endl;

	// Scan input for target
	array<int> o = search( fi, H, O.B);
	cout << argv[0] << ": searched." << endl;

	// Relabel silent parts
	if( O.thr > 0)
		for( size_t i = 0 ; i < o.size() ; ++i)
			if( sl(i) == 0)
				o(i) = -1;

	// Pass a median filter to keep sustained sections
	if( O.m){
#if 1
		cout << argv[0] << ": median-filtering state output." << endl;
		array<int> o2( o.size());
		for( size_t i = 0 ; i < o.size() ; ++i)
			o2(i) = o(i);
		int c[2];
		for( size_t i = O.m ; i < o.size()-O.m ; ++i){
			c[0] = c[1] = 0;
			for( int j = -O.m ; j <= O.m ; ++j)
				if( o2(i+j) != -1)
					++c[o2(i+j)];
			o(i) = O.mw*c[0] > c[1] ? 0 : 1;
		}
#else
		cout << argv[0] << ": lowpass-filtering state output." << endl;
		array<int> fo( o.size());
		for( int i = O.m ; i< o.size()-O.m ; ++i){
			fo(i) = 0;
			for( int j = -O.m ; j <= O.m ; ++j)
				fo(i) += o(i+j);
//			fo(i) += exp( -pow( 2.*j/O.m, 2)) * o(i+j);
			}
		for( int i = 0 ; i < o.size() ; ++i)
			o(i) = fo(i) > O.m/4;
#endif
	}

	// Show computation time
	gettimeofday( &time, NULL);
	const double end_time = time.tv_sec + double( time.tv_usec)/1000000;
	cout << argv[0] << ": features+learning+searching took " << end_time - start_time << " seconds." << endl;
	cout << argv[0] << ": processing was " << (double( sfi.frames)/sfi.samplerate) / (end_time - start_time) << " x real time." << endl;

	// Make an EDL list
	if( O.edl.size()){
		ofstream e( O.edl.c_str());
		bool cl = false;
		if( o(0) == 0){
			e << 0 << ' ';
			cl = true;
		}
		for( size_t i = 1 ; i < o.size() ; ++i)
			if( o(i) != o(i-1)){
				e << real_t(i*sz/O.hp)/sfi.samplerate;
				if( !cl){
					e << ' ';
					cl = true;
				}else{
					e << " 0\n";
					cl = false;
				}
			}
		if( cl)
			e << real_t( sfi.frames)/sfi.samplerate << " 0\n";
		cout << argv[0] << ": dumped EDL file " << O.edl << "." << endl;
	}

	// Make an output file with target sound
	int N = 0;
	wavfile_t sf( O.fo, ios::out, 1, sfi.samplerate);
	for( size_t i = 0 ; i < o.size() ; ++i)
		if( o(i) > stth){
			array<real_t> tt( in.v+i*int(sz/O.hp), sz/O.hp);
			sf.write( tt);
			++N;
		}
	cout << argv[0] << ": dumped " << N << " of " << o.size() << " frames to " << O.fo << "." << endl;
}
