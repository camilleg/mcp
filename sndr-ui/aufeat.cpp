#include "aufeat.h"
#include "optparse.h"
#include <iostream>

extern int ffmpeg_init();
extern bool ffmpeg_open(const std::string& filename);
extern bool ffmpeg_eof();
extern bool ffmpeg_read(array<double> &p);
extern double ffmpeg_samplerate();

#include "videostream.h"

int main(int argc, const char **argv)
{

	// Analysis options
	const std::string fopts = getoption<std::string>("-F", argc, argv, "cdm", "Feature options"); // e.g., a1w26c14dD
	const int fft_size = 1024;;;; // getoption<int>("-s", argc, argv, 1024, "FFT size"); // e.g., 512
	const int b = getoption<int>("-n", argc, argv, 13, "Feature coefficients"); // Number of coeffs
#ifdef UNUSED
	const real_t flo = getoption<real_t>("-l", argc, argv, 80, "Lowest frequency"); // in Hz
	const real_t fhi = getoption<real_t>("-h", argc, argv, 7500, "Highest frequency"); // in Hz
#endif
	const int fb = getoption<int>("-b", argc, argv, 32, "Filterbank filters");

	const std::string infileName = getoption<std::string>("-i", argc, argv, "/dev/null", "Input soundfile");
	const std::string outStream  = getoption<std::string>("-o", argc, argv, "/dev/null", "Output stream");

	if (ffmpeg_init() != 0)
	  return 1;

	// Open multimedia source (video or pure audio, filename or URL).
	if (!ffmpeg_open(infileName))
	  return 1;
	const real_t sr = ffmpeg_samplerate();

	// Allocate feature.
	aufeat_t<real_t> feature;
	feature.setup(sr, fft_size, fb, b, 50, sr/2, fopts);

	// Open video stream.
	videostream_init("127.0.0.1", "1234", outStream);

	array<real_t> x(fft_size);
	array<real_t> y(fft_size);
	array<real_t> Y;
	for (int i=0; !ffmpeg_eof(); ++i) {
		// Read buffer x from soundfile.
		if (!ffmpeg_read(x)) {
		  continue;
		}

		// Compute feature y from buffer x.
		feature.extract(y, x);

#if 0
		if (Y.empty()) {
#define hack 15000
			printf("resizing matrix Y to %lu by %d\n", y.size(), hack);
			//printf("2nd dim should be %d / #frames == i_max\n", fft_size); return 0;

			Y.resize(y.size(), hack); // fft_size/*wavfile.h's infile.frames*/ /fft_size+1;
		}
		printf("stuffing Y(%d .. %lu, %d)\n", 0, y.size()-1, i);
		for (size_t j = 0; j < y.size(); ++j)
			Y(j,i) = pow(y(j), 0.3);
#else
		// Export feature y to a video stream playable by ffplay.
		videostream_array(y);
#endif
	}
	std::cerr << std::endl;
	return 0;
}
