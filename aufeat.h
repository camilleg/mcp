// aufeat.h -- Class to extract audio features
// language: C++
// author  : Paris Smaragdis
// copyright: (c) 2013 University of Illinois, All Rights Reserved

#ifndef __AUFEAT_H__
#define __AUFEAT_H__

#include "array.h"
#include <fftw3.h> // apt-get install libfftw3-dev
#include <cmath>
#include <cfloat>
#include <fstream>

// Bark scale frequencies
const double Bark_freq[] = {
	50, 150, 250, 350, 450, 570, 700, 840, 1000, 
	1170, 1370, 1600, 1850, 2150, 2500, 2900, 3400, 
	4000, 4800, 5800, 7000, 8500, 10500, 13500
};

// Audio feature class
template <class T>
class aufeat_t{
public:
	T sr, min_f, max_f;
	size_t sz; int nf, b, o, od; // fft size, number of filters, number of coeffs, eventual output size, delta filter estimator size
	bool delta, delta2, zero, energy, bark, mel, dct, nrm, lg, tran;
	array<T> W, DCT, w, Tm, Tt;

	// Buffer crap due to delta estimation
	array<T> buf;
	int bc, cf;
	int cmod( int i, int j) const
	{
		while( i >= j)
			i -= j;
		while( i < 0)
			i += j;
		return i;
	}

	// FFT parameters
#ifdef __FLOAT
	fftwf_plan fft;
#else
	fftw_plan fft;
#endif
	array<T> fx;

	// Mel/Hz conversions
	T Hz2mel( T f) const { return 1127.01048 * log( 1.+f/700.); }
	T mel2Hz( T m) const { return 700.*(exp( m/1127.01048)-1.); }

	// Constructor/destructor
	aufeat_t() : sz(0) {}
	~aufeat_t()
	{
#ifdef __FLOAT
		fftwf_destroy_plan( fft);
#else
		fftw_destroy_plan( fft);
#endif
	}

	void report() const
	{
		std::cout << "pntr " << this << std::endl;
		std::cout << "param " << sr << ' ' << sz << ' ' << nf << ' ' << b << ' ' << min_f << ' ' << max_f << std::endl;
		std::cout << "buf " << buf.m << ' ' << buf.n << ' ' << buf.k << ' ' << buf.v << std::endl;
	}

	// Setup
	void setup( T _sr, int _sz, int _nf, int _b, T mn, T mx, const std::string &O, const std::string tf = "")
	{
		// Local copies
		sr = _sr;
		sz = _sz;
		nf = _nf;
		b = _b;
		min_f = mn;
		max_f = mx;

		// Parse the options string
		bark   = O.find( 'b') != std::string::npos;
		mel    = O.find( 'm') != std::string::npos;
		delta  = O.find( 'd') != std::string::npos;
		delta2 = O.find( 'D') != std::string::npos;
		zero   = O.find( '0') != std::string::npos;
		energy = O.find( 'e') != std::string::npos;
		dct    = O.find( 'c') != std::string::npos;
		nrm    = O.find( 'n') != std::string::npos;
		lg     = O.find( 'l') != std::string::npos;
		tran = tf.size() > 0;

		// Various checks
		if( mel && bark)
			std::cout << "Mel and bark at the same time?  I'll use bark then ..." << std::endl;
		if( max_f > sr/2){
			std::cout << "Maximum frequency (" << max_f << "Hz) more than sr/2 (" << sr/2 << "Hz), clipping to Nyquist instead" << std::endl;
			max_f = sr/2;
		}

		// Default output is the spectrogram so start with that output size
		o = sz/2+1;

		// Make a hann window
		w.resize( sz);
		for( size_t i = 0 ; i < sz ; ++i)
			w[i] = .5 + .5 * std::cos( -M_PI + 2.*M_PI*T(i)/sz);

		// Setup the FFT
		fx.resize( 2*(sz/2+2));
#ifdef __FLOAT
		fft = fftwf_plan_dft_r2c_1d( sz, &fx(0), (fftwf_complex*)&fx(0), FFTW_ESTIMATE);
#else
		fft = fftw_plan_dft_r2c_1d ( sz, &fx(0), (fftw_complex*)&fx(0), FFTW_ESTIMATE);
#endif

		// Center frequencies for warping
		array<T> fc;

		// Mel warping frequencies
		if( mel){
			if( min_f > 0){
				fc.resize( nf+1);
				T mn_f = Hz2mel( min_f);
				T mx_f = Hz2mel( max_f);
				for( int i = 0 ; i < nf ; i++)
					fc(i+1) = mel2Hz( mn_f + (mx_f-mn_f) * T(i)/(nf-1));
				fc(0) = 0;
			}else{
				fc.resize( nf);
				T mn_f = Hz2mel( min_f);
				T mx_f = Hz2mel( max_f);
				for( int i = 0 ; i < nf ; i++)
					fc[i] = mel2Hz( mn_f + (mx_f-mn_f) * T(i)/(nf-1));
			}
		}

		// Bark warping frequencies
		if( bark){
			// Find Barks inside our range
			int b1 = 0, b2 = 24;
			while( Bark_freq[b1] < min_f && b1 < 23) b1++;
			while( Bark_freq[b2-1] > max_f && b2-1 != 0) b2--;

			// Copy the list
			fc.resize( b2-b1+1);
			for( size_t i = 0 ; i < fc.size()-1 ; ++i)
				fc(i+1) = Bark_freq[b1+i];
			fc(0) = 0;
			nf = fc.size()-1;
		}

		// Make the warping matrix
		if( mel || bark){
			// Convert center frequencies to bin numbers
			for( size_t i = 0 ; i < fc.size() ; ++i)
				fc[i] = sz*fc[i]/sr;

			// Figure out the warping
			W.resize( sz/2+1, nf-1);
			for( size_t i = 1 ; i < fc.size()-1 ; ++i){
				T a1 = 1./(fc[i]-fc[i-1]), b1 =  fc[i-1]/(fc[i-1]-fc[i]);
				T a2 = 1./(fc[i]-fc[i+1]), b2 = -fc[i+1]/(fc[i]-fc[i+1]);
				for( size_t j = 0 ; j < sz/2+1 ; ++j){
					T as = a1*(j+1)+b1, ds = a2*(j+1)+b2; // why the j+1?  am i off by one? DC?
					as = (as > 1 || as < 0) ? 0 : as;
					ds = (ds > 1 || ds < 0) ? 0 : ds;
					W(j,i-1) = 2*(as + ds);
				}
			}
			nf--; // I'm skipping the last filter, what's its bandwidth?
			o = W.n;
//			std::ofstream ff( "warpdump.txt", std::ios::binary | std::ios::out);
//			ff << W << std::endl;
		}

		// Load the transformation data
		if( tran){
			std::ifstream f( tf.c_str(), std::ios::in | std::ios::binary);
			int wm_wn[2];
			f.read( (char*)wm_wn, 2 * sizeof( int));
			Tt.resize( wm_wn[0], wm_wn[1]);
			f.read( (char*)Tt.v, Tt.size()*sizeof( T));
			Tm.resize( wm_wn[0]);
			f.read( (char*)Tm.v, Tm.size()*sizeof( T));
			o = Tt.n;
		}

		// Make the DCT transform
		if( dct){
			DCT.resize( b-!zero, nf);
			for( int i = !zero ; i < b ; i++)
				for( int j = 0 ; j < nf ; j++)
					DCT(i-!zero,j) = std::sqrt( 2./nf) * std::cos( (M_PI/nf) * (j+.5) * i);
			if( zero)
				for( int j = 0 ; j < nf ; j++)
					DCT(0,j) /= std::sqrt( 2.); // that's what matlab returns for a dct ...
			o = DCT.m;
//			std::ofstream ff( "dctdump.txt", std::ios::binary | std::ios::out);
//			ff << DCT << std::endl;
		}

		// Make the delta filter flags and constants
		if( delta){
			od = 4;
			o *= 2;
		}else
			od = 0;

		// Setup the buffer for the delay filter and init it
//		cout << buf.m << ' ' << buf.n << ' ' << buf.k << ' ' << buf.v << endl;;
//		cout << o << ' ' << 2*od+1 << endl;
		buf.resize( o, 2*od+1);
//		cout << o << ' ' << 2*od+1 << endl;
		for( size_t i = 0 ; i < buf.size() ; ++i)
			buf[i] = 0;
		bc = 0;
		cf = od;
	}

	// Extract features from an input, if using delta features it imposes a delay of 4 frames
	T extract( array<T> &y, const array<T> &x, const bool use_buf = true)
	{
		using namespace std;
#ifdef __CHECK
		// Check sizes
		if( x.size() != sz)
			throw std::runtime_error( "audfeat_t::extract(): Input array's size differs from transform size");
#endif
		// Setup and do the FFT
		for( size_t i = 0 ; i < sz ; ++i)
			fx[i] = w[i] * x[i] + 1e-6;
#ifdef __FLOAT
		fftwf_execute( fft);
#else
		fftw_execute( fft);
#endif
		// Take the magnitude
		for( size_t i = 0 ; i < sz/2+1 ; ++i)
			fx[i] = sqrt( fx(2*i)*fx(2*i) + fx(2*i+1)*fx(2*i+1));

		if( nrm){
			// Normalize
			T sm = FLT_EPSILON;
			for( size_t i = 0 ; i < sz/2+1 ; ++i)
				sm += fx[i];
			for( size_t i = 0 ; i < sz/2+1 ; ++i)
				fx[i] /= sm;
		}

		if( lg){
			// Take the log
			for( size_t i = 0 ; i < sz/2+1 ; ++i)
				fx[i] = log( fx[i] + FLT_EPSILON);
		}

		// Warp according to mapping
		array<T> t;
		if( mel || bark){
			t.resize( W.n);
			for( size_t i = 0 ; i < W.n ; ++i){
				T tt = 0;
				for( size_t j = 0 ; j < W.m ; ++j)
					tt += W(j,i)*fx(j);
				t[i] = log( tt + FLT_EPSILON);
			}
		}else if( tran){
			t.resize( Tt.n);
			for( size_t i = 0 ; i < Tt.n ; ++i){
				T tt = 0;
				for( size_t j = 0 ; j < Tt.m ; ++j)
					tt += Tt(j,i)*(fx(j)-Tm(j));
				t[i] = tt;
			}
		}else{
			t.resize( sz/2+1);
			for( size_t i = 0 ; i < t.size() ; ++i)
				t[i] = fx[i];
		}

		// DCT it
		array<T> t2;
		if( dct){
			t2.resize( DCT.m);
			for( size_t i = 0 ; i < DCT.m ; ++i){
				T tt = 0;
				for( size_t j = 0 ; j < DCT.n ; ++j)
					tt += DCT(i,j)*t(j);
				t2[i] = tt;
			}
		}else{
			t2.resize( t.size());
			for( size_t i = 0 ; i < t.size() ; ++i)
				t2[i] = t[i];
		}

		// Get input window's energy
		T sm = 0;
		for( size_t i = 0 ; i < x.size() ; ++i)
			sm += fabs( x[i]);
		sm /= x.size();

		// Hack to avoid NaNs/Infs
		if( sm == 0)
			for( size_t i = 0 ; i < t2.size() ; ++i)
				t2[i] = 0;

		if( use_buf){
			// Add to buffer
			for( size_t i = 0 ; i < t2.size() ; ++i){
				if( delta){
					buf(t2.size()+i,cf) = 0;

					// Filter to get deltas
					for( int j = -od ; j <= od ; j++)
						buf(t2.size()+i,cf) += j*buf(i,cmod(cf+j,buf.n))/60;
				}
				buf(i,bc) = t2[i];
			}

			// Copy result to output
			y.resize( o);
			for( size_t i = 0 ; i < buf.m ; ++i)
				y[i] = buf(i,cf);

			// Update buffer counters
			++bc %= buf.n;
			++cf %= buf.n;
		}else{
			// Don't use buffer (or deltas)
			y.resize( t2.size());
			for( size_t i = 0 ; i < y.size() ; ++i)
				y[i] = t2[i];
		}

		return sm; // frame energy
	}

	// Extract features from an input, offline with no delta delay
	void extract_offline( array<T> &y, array<T> &e, const array<T> &x, const int hp)
	{
		y.resize( o, (x.size()-sz)/hp + 1);
		e.resize(    (x.size()-sz)/hp + 1);
		int ri = 0;
		size_t to = 0;
		for( int i = 0 ; i < int(x.size()-sz) ; i+=hp,++ri){
			// Make temps
			array<T> tx( x.v+i, sz), ty; //( y.v+ri*o, o);

			// Get non-delta features
			e(ri) = extract( ty, tx, false);
			to = ty.size();
			for( size_t j = 0 ; j < to ; ++j)
				y(j,ri) = ty(j);
		}

		if( delta)
			// Compute the deltas
			for( size_t i = 0 ; i < y.n ; ++i)
				for( size_t k = 0 ; k < to ; ++k){
					y(to+k,i) = 0;
					for( int j = -od ; j <= od ; ++j)
						y(to+k,i) += j * y(k,std::min( std::max( size_t(0), i+j), y.n-1)) / 60;
				}
	}
};

#endif
