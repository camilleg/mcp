// sndr.cpp -- Driver source code for sound recognition tasks
// language: C++
// author  : Paris Smaragdis

#include "sndr.h"

// Timer
double _tic_start_time;
// Start counting
void tic()
{
	struct timeval time;
	gettimeofday( &time, NULL);
	_tic_start_time = time.tv_sec + double( time.tv_usec)/1000000.0;
}
// Get elapsed time
double toc()
{
	struct timeval time;
	gettimeofday( &time, NULL);
	const double end_time = time.tv_sec + double( time.tv_usec)/1000000.0;
	return end_time - _tic_start_time;
}

int main( int argc, const char **argv)
{
	// Get options

	// Filenames
	std::string modout = getoption<std::string>( "-M", argc, argv, "", "Output model filename");
	const std::string dump = getoption<std::string>( "-d", argc, argv, "", "Dump file prefix");
	const array<std::string> infile = mgetoption<std::string>( "-i", argc, argv, "Input soundfile(s)");
	const array<std::string> modin = mgetoption<std::string>( "-m", argc, argv, "Input model(s)");
	const std::string edl = getoption<std::string>( "-D", argc, argv, "", "EDL filename prefix");
	const array<std::string> target = mgetoption<std::string>( "-g", argc, argv, "Target file(s)"); // enables sound tracking

	// Model
	const int K = getoption<int>( "-K", argc, argv, 8, "Number of Gaussians"); // per GMM
	const int it = getoption<int>( "-e", argc, argv, 50, "Number of learning iterations");
	const real_t trans = getoption<real_t>( "-p", argc, argv, -40, "Transition likelihood"); // log lik for other states if negative, lik for diagonal if positive
	const array<real_t> r = mgetoption<real_t>( "-r", argc, argv, "State transition bias");
	const array<real_t> bias = mgetoption<real_t>( "-B", argc, argv, "Likelihood bias");
	const int mf = getoption<int>( "-f", argc, argv, 0, "Length of state output filter");
	const real_t mw = getoption<real_t>( "-w", argc, argv, 1, "Bias of state output filter");

	// Analysis
	AudioFeatureExtractor_t<real_t> A;
	A.fopts = getoption<std::string>( "-F", argc, argv, "cdm", "Feature options");
	A.b = getoption<int>( "-n", argc, argv, 13, "Feature coefficients"); // Number of coeffs to use
	A.flo = getoption<real_t>( "-l", argc, argv, 80, "Lowest frequency");
	A.fhi = getoption<real_t>( "-h", argc, argv, 7500, "Highest frequency");
	A.fb = getoption<int>( "-b", argc, argv, 32, "Filterbank filters");

	// Time and amplitude scale
	A.thr = getoption<real_t>( "-T", argc, argv, 0, "Peak threshold");
	A.tsz = getoption<real_t>( "-t", argc, argv, 0.1, "Window size"); // in seconds
	A.hp = getoption<int>( "-H", argc, argv, 1, "Hop size"); // window size/hp
	A.av = getoption<int>( "-a", argc, argv, 1, "Feature averaging");

	tic();

	// Choose whether to train, classify, or track.
	if( !modin.empty()){

		// Classify an input sound from pre-trained models
		std::cout << argv[0] << ": classifying..." << std::endl;
		AudioClassifier_t<real_t> C( A);

		// Set some decoding parameters
		C.filter( mf, mw);
		C.trans = trans;
		C.bias = bias;
		C.r = r;

		// Load classifiers and combine them if necessary
		if( modin.size() == 1) {
			C.load( modin(0)); // may throw exception
		} else {
			if (modout.empty()) {
				modout = "model-default";
				std::cout << argv[0] << ": output model filename defaulting to " << modout << "." << std::endl;
			}
			C.combine_models( modin, modout); // may throw exception
		}

		// For each file...
		for( size_t i = 0 ; i < infile.size() ; ++i){
			wavfile_t f( infile(i), std::ios::in);
			array<real_t> x( f.frames);
			f.read_mono( x);

			// ...perform classification...
			C( x, int( f.samplerate));

			// ... and output .wav and .edl files.
			if (dump != "") {
				for (int j= 0; j<C.H.S; ++j)
					C.make_snd( x, infile(i) + "." + dump + "." + to_str(j) + ".wav", j);
			}
			if (edl != "") {
				for (int j=0; j<C.H.S; ++j)
					C.make_edl( infile(i) + "." + edl + "." + to_str(j) + ".edl", j);
			}

		}

	}else if( !target.empty()){

		// Track a sound
		std::cout << argv[0] << ": learning..." << std::endl;
		
		// Add example sounds from the dictionary
		AudioFeatures_t<real_t> F( A);
		for( size_t i = 0 ; i < infile.size() ; ++i){
			std::cout << argv[0] << ": learning example sound " << infile(i) << std::endl;
			wavfile_t f( infile(i), std::ios::in);
			array<real_t> x( f.frames);
			f.read_mono( x);
			F( x, f.samplerate, true);
		}
		
		// Learn and save model
		AudioModel_t<real_t> M( A, K);
		M( F, it);
		M.save( modout + "-ubm");

		// Add target sounds from dictionary
		F.clear();
		for( size_t i = 0 ; i < target.size() ; ++i){
			std::cout << argv[0] << ": learning target sound " << target(i) << std::endl;
			wavfile_t f( target(i), std::ios::in);
			array<real_t> x( f.frames);
			f.read_mono( x);
			F( x, f.samplerate, true);
		}

		// Learn and save target file based on the universal model
		AudioModel_t<real_t> Mt( A, K);
		Mt( F, it, M.G);
		Mt.save( modout + "-target");

	}else{
		// Learn a model from example sounds
		std::cout << argv[0] << ": learning..." << std::endl;

		// Add example sounds from dictionary
		AudioFeatures_t<real_t> F( A);
		for( size_t i = 0 ; i < infile.size() ; ++i){
			wavfile_t f( infile(i), std::ios::in);
			array<real_t> x( f.frames);
			f.read_mono( x);
			F( x, f.samplerate, true);
		}

		// Learn and save model
		AudioModel_t<real_t> M( A, K);
		M( F, it);
		M.save( modout);
	}

	std::cout << argv[0] << ": done in " << toc() << " seconds." << std::endl;
	return 0;
}
