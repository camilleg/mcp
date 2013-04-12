#include "wavfile.h"
#include "aufeat.h"
#include "optparse.h"

#include <iostream>

#include "aquat.h"

using namespace std;

typedef double real_t;

int main( int argc, const char **argv)
{
	// Analysis options
	const string fopts = getoption<string>( "-F", argc, argv, string( "cdm"), "Feature options"); // e.g., a1w26c14dD
	const int sz = getoption<int>( "-s", argc, argv, 1024, "FFT size"); // e.g., 512
	const int b = getoption<int>( "-n", argc, argv, 13, "Feature coefficients"); // Number of coeffs
#ifdef UNUSED
	const real_t flo = getoption<real_t>( "-l", argc, argv, 80, "Lowest frequency"); // in Hz
	const real_t fhi = getoption<real_t>( "-h", argc, argv, 7500, "Highest frequency"); // in Hz
#endif
	const int fb = getoption<int>( "-b", argc, argv, 32, "Filterbank filters");

	const string infile = getoption<string>( "-i", argc, argv, string( "80s.wav"), "Input soundfile");

	// Open the soundfile
	wavfile_t w( infile);
	const real_t sr = w.samplerate;

	// Allocate features object
	aufeat_t<real_t> f;
	f.setup( sr, sz, fb, b, 50, sr/2, fopts);

	// Open a graphics terminal
	aq_window( 1, "Features");

	// Process
	array<real_t> x( 1024);
	array<real_t> y( 1024);
	array<real_t> Y;
	for(int i=0; w.file; ++i){
		cout << "." << endl;
		// Read a buffer
		w.read_mono( x);

		// Get its features
		f.extract( y, x);

		// Store them inside a matrix
		if( !Y.size())
			Y.resize( y.size(), w.frames/sz+1);
		for( size_t j = 0 ; j < y.size() ; ++j)
			Y(j,i) = pow( y(j), 0.3);
	}

	// Show me
	aq_image( Y.v, Y.m, Y.n);
}
