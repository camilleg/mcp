// sndr.h -- Classes to perform sound classification tasks
// language: C++
// author  : Paris Smaragdis

#ifndef __SNDR_H
#define __SNDR_H

#include "intransp.h"
#include "aufeat.h"
#include "hmm.h"
#include "optparse.h"
#include "wavfile.h"

#include <list>
#include <sys/time.h>
#include <vector>

template <class T>
class AudioFeatureExtractor_t {
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

  // Constructors
  AudioFeatureExtractor_t(int _b=13, T _tsz=.1, T _flo=120, T _fhi=6855,
      int _fb=30, std::string _fopts=std::string("cdm"), T _thr=0, int _av=1) : 
    b(_b), fb(_fb), hp(1), av(_av), sz(0), tsz(_tsz),
    flo(_flo), fhi(_fhi), thr(_thr), srate(0), fopts(_fopts) {}

  AudioFeatureExtractor_t(AudioFeatureExtractor_t &f) : 
    b(f.b), fb(f.fb), hp(f.hp), av(f.av), sz(0), tsz(f.tsz),
    flo(f.flo), fhi(f.fhi), thr(f.thr), srate(0), fopts(f.fopts), F(f.F) {}

  void report() {
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

  // Sneakily append feature data to another file.
  void save(const std::string& filename) {
    std::ofstream ff(filename.c_str(), std::ios::binary | std::ios::app);
    ff.seekp(0, std::ios::end);
    ff.write("feat", 4*sizeof(char));
    ff.write((char*)&b, sizeof(int));
    ff.write((char*)&fb, sizeof(int));
    ff.write((char*)&hp, sizeof(int));
    ff.write((char*)&av, sizeof(int));
    ff.write((char*)&tsz, sizeof(T));
    ff.write((char*)&flo, sizeof(T));
    ff.write((char*)&fhi, sizeof(T));
    ff.write((char*)&thr, sizeof(T));
    ff.write((char*)&srate, sizeof(T));
//  // Member array<T> bias was moved to class AudioClassifier_t.
//  const int bs = bias.size();
//  ff.write((char*)&bs, sizeof(int));
//  ff.write((char*)bias.v, bias.size()*sizeof(T));
    ff.write(fopts.c_str(), fopts.size()*sizeof(char));
  }

  // Search for sneakily appended feature data.
  void load(const std::string& filename) {
    char opt[32];
    char *opp = opt;
    std::ifstream ff(filename.c_str(), std::ios::binary);
    do {
      ff.read(opt, 4*sizeof(char));
      ff.seekg(-3, std::ios::cur);
    } while (!ff.eof() && !(opt[0]=='f' && opt[1]=='e' && opt[2]=='a' && opt[3]=='t'));
    opt[4] = '\0';
    if (!strcmp(opt, "feat")){
      std::cout << "Found feature parameters appended to file " << filename << "." << std::endl;
      ff.seekg(3, std::ios::cur);
      ff.read((char*)&b, sizeof(int));
      ff.read((char*)&fb, sizeof(int));
      ff.read((char*)&hp, sizeof(int));
      ff.read((char*)&av, sizeof(int));
      ff.read((char*)&tsz, sizeof(T));
      ff.read((char*)&flo, sizeof(T));
      ff.read((char*)&fhi, sizeof(T));
      ff.read((char*)&thr, sizeof(T));
      ff.read((char*)&srate, sizeof(T));
//    // Member array<T> bias was moved to class AudioClassifier_t.
//    int tbs;
//    ff.read((char*)&tbs, sizeof(int));
//    if (tbs > 0){
//	bias.resize(tbs);
//	ff.read((char*)bias.v, bias.size()*sizeof(T));
//    }
      while (!ff.eof()) ff.read(opp++, sizeof(char));
      *--opp = '\0';
      fopts = opt;
    }
  }

  // Check options match
  bool operator==(const AudioFeatureExtractor_t<T> &A) const
  {
    return b == A.b && fb == A.fb && hp == A.hp && av == A.av &&
	    tsz == A.tsz && flo == A.flo && fhi == A.fhi && thr == A.thr &&
	    srate == A.srate && fopts.compare(A.fopts) == 0;
  }
  bool operator!=(const AudioFeatureExtractor_t<T> &A) const { return !(*this == A); }
								  
  // Extract features f and peaks p from sound s.
  void operator()(array<T> &f, array<int> &p, const array<T> &s, const int sr, const bool thrm=false)
  {
    // Init the feature structure according to the new sample rate
    if (sr != srate){
      sz = pow(2., int(log2(tsz*sr)));
      F.setup(sr, sz, round(3*log(sr)), b, flo, fhi, fopts);
//    F.report();
      srate = sr;
    }

    // Some preallocation
    array<T> e((s.size()-sz)/(sz/hp)+1);
    
#if 0
    // Get features of input sound using online estimation (has 4-frame delay)
    f.resize (F.o, (s.size()-sz)/(sz/hp)+1);
    for (int i=0,j=0; i<s.size()-sz; i+=sz/hp, j++){
      array<T> st(s.v+i, sz), ft(f.v+j*F.o, F.o);
      e(j) = F.extract(st, ft);
    }
#else
    // Get features of input sound using offline estimation
    F.extract_offline(f, e, s, sz/hp);
#endif

    if (thr > 0.0) {
      // Threshold is nonzero, so remove low-energy frames.
      const T peakVolume = thr * e.max();
      
      // Mark non-quiet frames
      p.resize(e.size());
      size_t cPeak = 0;
      for (size_t i=0; i<e.size(); ++i) {
	p(i) = e(i) >= peakVolume;
	if (p(i)) ++cPeak;
      }

      if (thrm) {
	// Remove quiet frames.
	// TODO: replace this copying with an in-place copy.  Might be faster.
	const array<T> f2(f);

	// Keep only the loud parts
	f.resize(F.o, cPeak);
	for (size_t i=0,j=0; i<e.size(); ++i) {
	  if (p(i)) {
	    for (int k=0; k<F.o; ++k) f(k,j) = f2(k,i);
	    ++j;
	  }
	}
	std::cout << "Volume trimmed from " << f2.n << " to " << f.n << " frames." << std::endl;
      }
}
    
    if (av > 1) {
      // Feature averaging
      const array<T> f2(f);
      f.resize(F.o, f.n/av);
      for (size_t i=0; i<f.m; ++i)
	for (size_t j=0; j<f.n; ++j) {
	  T t = 0;
	  for (int k=0; k<av; ++k)
	    if (j+k < f2.n)
	      t += f2(i, j+k);
	  f(i,j) = t/av;
	}
    }

    // Transpose in place to make the cache happy during training
    // (How much does that matter, when f is entirely copied a moment before?)
    intp(&f[0], f.m, f.n);
    f.k=f.m; f.m=f.n; f.n=f.k; f.k=1;

#if 0
    // Dump to debugging file
    {
      std::ofstream df("/tmp/sndr-debugging-dump.dat", std::ios::binary | std::ios::out);
      df.write((char*)&f.m, sizeof(int));
      df.write((char*)&f.n, sizeof(int));
      df.write((char*)f.v, sizeof(T)*f.m*f.n);
    }
#endif
  }
};

// Audio feature container
template <class T>
class AudioFeatures_t {
public:
  AudioFeatureExtractor_t<T> &F; // Feature extractor to use
  array<T> C; // Consolidated features
  std::list<array<T> > D; // Feature data
  std::list<array<int> > S; // Silent frame data

  // Initialize from a feature extractor
  AudioFeatures_t(AudioFeatureExtractor_t<T> &_F) : F(_F) {}

  // Append sound "in" to learning dictionary D and S
  void operator()(const array<T> &in, const T sr, bool thrm=false)
  {
    D.push_back(array<T>());
    S.push_back(array<int>());
    F(D.back(), S.back(), in, sr, thrm);
  }

  // Consolidate all feature sets
  void consolidate()
  {
    // Count our data
    int fs = 0;
    for (typename std::list<array<T> >::iterator i = D.begin(); i != D.end(); ++i)
      fs += (*i).m;

//  // Nothing to do
//  if (fs == 0)
//    throw std::runtime_error("AudioFeatures_t<T>::consolidate(): empty list");

    // Consolidate features
    int ck = 0;
    const int del = F.fopts.find('d') != std::string::npos;
    C.resize(fs-del*5*D.size(), D.front().n);
    while (D.size()) {
      for (unsigned int k = del*5; k < D.front().m; ++k, ++ck)
	for (unsigned int j = 0; j < D.front().n; ++j)
	  C(ck,j) = D.front()(k,j);
      D.pop_front();
      S.pop_front();
    }
    std::cout << "Overall consolidated feature size is " << C.m << " x " << C.n << "." << std::endl;
  }

  // Clear buffers
  void clear()
  {
    C.resize(0);
    while (!D.empty()) D.pop_front();
    while (!S.empty()) S.pop_front();
  }
};

// Single-class model
template <class T>
class AudioModel_t {
public:
#ifdef __HMM_TRAIN
  #define baseModel_t hmm_t
#else
  #define baseModel_t gmm_t
#endif
  baseModel_t G;                 // HMM or GMM model of sound class
  AudioFeatureExtractor_t<T> &F; // Feature extractor reference, used as a test

  // Constructors
  AudioModel_t(AudioFeatureExtractor_t<T> &_F): F(_F) {}
  AudioModel_t(AudioFeatureExtractor_t<T> &_F, const int k): G(k), F(_F) {}

  // Train model G from data A.
  void operator()(AudioFeatures_t<T> &A, const int it) {
    if (A.F != F)
      throw std::runtime_error("AudioModel_t trained feature space mismatch.");
    A.consolidate();
    G.train(A.C, it);
    std::cout << "Trained model in " << it << " iterations." << std::endl;
  }
  
  // Train model G from data A and an initial model.
  void operator()(AudioFeatures_t<T> &A, const int it, const baseModel_t Ginitial)
  {
    if (A.F != F)
      throw std::runtime_error("AudioModel_t model-trained feature space mismatch.");
    A.consolidate();
    G.train(A.C, it, Ginitial);
    std::cout << "Updated model in " << it << " iterations." << std::endl;
  }
  
  bool save(const std::string &f) const {
    if (f.empty()) {
      std::cout << "error: AudioModel_t::save(\"\");" << std::endl;
      return false;
    }
    G.save(f);
    F.save(f); // Append feature parameters.
    std::cout << "Saved AudioModel_t model " << f << "." << std::endl;
    return true;
  }

  bool load(const std::string &f) {
    if (f.empty()) {
      std::cout << "error: AudioModel_t::load(\"\");" << std::endl;
      return false;
    }
    G.load(f);
    F.load(f); // Load any appended feature parameters.

    F.srate = 0; // Ensure feature class initialization.
    std::cout << "Loaded model " << f << "." << std::endl;
    return true;
  }
};

// Multiple-class combiner
template <class T>
class AudioClassifier_t {
public:
  int m;  // Length of state output median filter
  T mw;   // Median filter bias
  
  T trans;       // HMM transition likelihood	
  array<T> r;    // Transition bias between states
  array<T> bias; // Likelihood bias for each HMM state

  array<int> o;  // Classification output
  hmm_t H;       // Master HMM
  AudioFeatureExtractor_t<T> &F;

  // Initialize
  AudioClassifier_t(AudioFeatureExtractor_t<T> &_F) : m(0), mw(1.), trans(-60), F(_F) {}

  // Combine multiple saved sound models into one classifier
  void combine(const array<std::string > &f) {
    std::list<AudioModel_t<T> > Al; // accumulator for models
    for (size_t i=0; i<f.size(); ++i) {
      if (f(i).empty()) {
	std::cout << "AudioClassifier_t::combine() skipping empty filename." << std::endl;
	continue;
      }
      Al.push_back(AudioModel_t<T>(F));
      if (!Al.back().load(f(i))) {
	std::cout << "AudioClassifier_t::combine() failed to load model '" << f(i) << "'." << std::endl;
      }
    }
    combine(Al);
  }

  // Combine multiple sound models into one classifier
  void combine(const std::list<AudioModel_t<T> > &Al)
  {
#ifdef __HMM_TRAIN
    // Pack into a single HMM.
    if (Al.size() != 2)
      throw std::runtime_error("BUG: hmm_t::combine() ignores intermediate elements of list of HMMs.");
    // TODO: use all elements of list, not just the first and last elements.
    ::combine(H, Al.front().G, Al.back().G, T(1.-trans), T(1.-trans));
#else
    // Pack all GMMs in a HMM
    int ai = 0;
    for (typename std::list<AudioModel_t<T> >::const_iterator A=Al.begin(); A!=Al.end(); ++A,++ai) {
//    // Make sure all models are relevant
//    if ((*A).F != F)
//      throw std::runtime_error("AudioClassifier_t<T>::combine(): Input models are not using the same features");

      // Get the model's GMM
      const gmm_t& G = A->G;

      if (ai == 0) {
	// First time through loop.  Allocate memory in HMM.
	H.S = Al.size();
	H.K = G.numGaussians();
	std::cout << "Making a " << H.S << "-state HMM with " << H.K << " gaussians per state." << std::endl;
	H.lPi.resize(H.S);
	H.lA.resize(H.S, H.S);
	H.smps.resize(H.S);
      }

      // Copy GMMs into HMM.
      H.smps(ai).copyGuts(G);

      // Make the transition matrix row
      if (trans > 0){
	for (int i=0; i<H.S; ++i) H.lA(ai,i) = log((1.-trans)/(H.S-1));
	H.lA(ai,ai) = log(trans);
      }else{
	for (int i=0; i<H.S; ++i) H.lA(ai,i) = trans;
	H.lA(ai,ai) = log(0.9999); // **** improper value, should be ok further down when I normalize lA
      }			
  }

    // Set initial probabilities H.lPi.
    for (int s=0; s<H.S; ++s) H.lPi(s) = log(1.0/H.S);
    
    for (int s=0; s<H.S; ++s) H.smps(s).normalizeLogPriors();

    // Bias the transitions H.lA.
    const int c = H.S * H.S;
    if (int(r.size()) == c) {
      std::cout << "Biasing" << std::endl;
      for (int i=0; i<c; ++i) H.lA(i) -= r(i);
    }

    // Re-normalize the transition matrix H.lA.
    for (int i=0; i<H.S; ++i){
      T ls = log(0.0);
      for (int j=0; j<H.S; ++j) ls = std::max(ls, H.lA(i,j)) + log1p(exp(-fabs(ls - H.lA(i,j))));
      for (int j=0; j<H.S; ++j) H.lA(i,j) -= ls;
    }
#endif
    std::cout << "Transition matrix:\n" << H.lA;
  }

  // Classify a new input
  void operator()(const array<T> &in, const int sr)
  {
    // Get the sound features
    array<T> D;
    array<int> S;
    F(D, S, in, sr, false);

    // Classify
    H.classify(o, D, bias);
    std::cout << "Input length " << in.size() << ", window size " << F.sz << ", output length " << o.size() << "." << std::endl;

    // Relabel the silent parts to the last known class
    if (F.thr > 0) {
      int c = 0;
      for (size_t i=0; i<o.size(); ++i) {
	if (S(i))
	  c = o(i);
	else
	  o(i) = c;
      }
    }

    // Pass through a median filter to smooth out sustained sections
    if (m > 1) {
      array<int> o2(o);
      // Don't use size_t i, which can wrap around to "negative" values.
      for (int i=m; i < int(o.size()-m); ++i) {
	if (i < 0) {
	  // Maybe m is preposterously large, like larger than o.size().
	  std::cout << "Median filter tried to misbehave: " << i << "," << o.size() << "," << m << "." << std::endl;
	  break;
	}
	int c[2] = {0};
	for (int j = -m; j <= m; ++j)
	  if (o2(i+j) != -1)
	    ++c[o2(i+j)];
	o(i) = mw*c[0] > c[1] ? 0 : 1;
      }
    }

    // Show the rates
    array<int> hs(H.S);
    for (int j=0; j != H.S; ++j)
      hs(j) = 0;
    for (size_t j=0; j != o.size(); ++j)
      hs(o(j))++;
    std::cout << "Class result histogram: " << hs << std::endl;
  }

  // Set state filtering parameters
  void filter(const int _m, const T _w)
  {
    m = _m;
    mw = _w;
  }

  // Make an EDL file for class ii
  void make_edl(const std::string &edl, const int ii=0) {
    std::ofstream e(edl.c_str());
    bool cl = false;
    if (o(0) == ii) {
      e << 0 << ' ';
      cl = true;
    }
    for (size_t i=1; i < o.size(); ++i) {
      if (o(i) != o(i-1)) {
	e << T(i*F.sz/F.hp) / F.srate;
	if (!cl && o(i) == ii) {
	  e << ' ';
	  cl = true;
	} else {
	  e << " 0\n";
	  cl = false;
	}
      }
    }
    if (cl)
      e << T(o.size()*F.sz/F.hp) / F.srate << " 0\n";
  }

  // Make a soundfile for class ii
  void make_snd(const array<T> &x, const std::string &f, const int ii=0) {
    wavfile_t sf(f, std::ios::out, 1, F.srate);
    int N=0;
    for (size_t k=0; k<o.size(); ++k) {
      if (o(k) == ii) {
	sf.write(array<T>(x.v+k*F.sz/F.hp, F.sz/F.hp));
	++N;
      }
    }
    std::cout << "Dumped " << N << " of " << o.size() << " frames to " << f << "." << std::endl;
  }

  // Load HMM model
  void load(const std::string& filename) {
    // Load the HMM, and any appended feature parameters.  Throw exception on failure.
    H.load(filename);
    F.load(filename);

    // Ensure feature class initialization.
    F.srate = 0;
  }

  // Combine multiple sound model files (and feature info) into an HMM model file.
  void combine_models(const array<std::string>& filenames_src, const std::string& filename_dst) {
    // Load GMM models and combine them into one HMM
    combine(filenames_src);

    // Save the HMM.  Append feature parameters.  Throw exception on failure.
    H.save(filename_dst);
    F.save(filename_dst);
  }
};

#endif
