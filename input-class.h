#include <cstdlib> // for rand()
#include <cassert>

class stream_t {
};

class input_t {
public:

  enum sample_format_t {
    // Copied from /usr/include/sndfile.h, from apt-get install libsndfile1-dev.
    SF_FORMAT_WAV                   = 0x010000,
    SF_FORMAT_AIFF                  = 0x020000,
    SF_FORMAT_PCM_16                = 0x0002,
    SF_FORMAT_ULAW                  = 0x0010,
    SF_ENDIAN_LITTLE                = 0x10000000
  };

  // Characteristics, aka attributes, aka specification, of a stream.
  class specifications_t {
  public:
    // These are const in spirit but not in fact,
    // merely so the implicit assignment op works.
    int channels;
    double sample_rate;
    enum sample_format_t sample_format;

    specifications_t(int c=1, double r=48000.0, sample_format_t f=SF_FORMAT_WAV) :
      channels(c), sample_rate(r), sample_format(f)
      { assert(sample_rate > 0.0); }
  };
  specifications_t specifications;

  // Return a fresh chunk of signal, as channels by sample frames.
  // (Because size_t can convert to double, overloading one read() with the two obvious signatures would be ambiguous.)
  // channel_mask is a bitmask, e.g. 0b00111111 for the first 6 channels.
  // (Later, eigen<> might work as easily as array<>.)
  array<short>& readSamples( size_t n, size_t offset, int channel_mask);
  array<short>& readSeconds( double n, double offset, int channel_mask);

  // Constructors and destructors.
  input_t() {}
  ~input_t() {} // close file, socket, etc.
  // Copy stream's specification
  input_t(const input_t& rhs, stream_t stream) :
    specifications(rhs.specifications) {}
  input_t(const input_t::specifications_t rhs, stream_t stream) :
    specifications(rhs) {}

  // Utilities
  double sample_rate() const {
    const double& r = specifications.sample_rate;
    assert(r > 0.0);
    return r;
  }
  int channels() const { return specifications.channels; }
  enum sample_format_t sample_format() const { return specifications.sample_format; }
  bool eof() const { return false; }
  operator bool() const { return true; } // still valid
};

array<short>& input_t::readSamples(size_t n, size_t offset_samples, int channel_mask) {
  if ((long)n < 0) {
    std::cerr << "input_t::readSamples zeroing negative sample count " << (long)n << ".\n";
    n = 0;
  }
  // If n==0, return an empty array.
  array<short>* a = new array<short>(n);
  for (size_t i=0; i<n; ++i)
    (*a)(i) = short(drand48()*65535 - 32768);
  return *a;
}

array<short>& input_t::readSeconds(double seconds, double offset_sec, int channel_mask) {
  if (seconds < 0.0) {
    std::cerr << "input_t::readSeconds zeroing negative seconds-duration " << seconds << ".\n";
    seconds = 0.0;
  }
  const int n = seconds * sample_rate();
  std::cerr << n << " samples\n";
  return readSamples(seconds * sample_rate(), offset_sec * sample_rate(), channel_mask);
}
