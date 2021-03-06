// wavfile.h -- Read/write for WAV files
// language: C++
// author  : Paris Smaragdis
// copyright: (c) 2009 Adobe Systems Incorporated, All Rights Reserved

#pragma once

#include "array.h"
#include <fstream>
#include <string>

class wavfile_t{
public:
	unsigned int channels;
	unsigned int frames;
	double samplerate;
	std::fstream file;
	int file_type;
	int bits;

	// Open a wave file
	wavfile_t( const std::string &filename, std::ios::openmode mode = std::ios::in, unsigned int ch = 1, double sr = 16000) 
	 : channels( ch), frames( 0), samplerate( sr), file( NULL), file_type( 0), bits( 0)
	{
		if (filename.empty())
			throw std::runtime_error( "wavfile_t::wavfile_t(): empty filename.");

		// Open an input file
		if( mode == std::ios::in){
			file.open( filename.c_str(), std::ios::in | std::ios::binary);
			if (!file)
				throw std::runtime_error( "wavfile_t::wavfile_t(): failed to read '" + filename + "'.");

			// at start of file
			file.seekg( 0, std::ios::beg);

			// Check if file is in WAVE format
			char id[4];
			file.read( id, 4);
			if( memcmp( id, "RIFF", 4))
				throw std::runtime_error( "wavfile_t::wavfile_t(): no starting 'RIFF' in file '" + filename + "'.");

			file.seekg( 4, std::ios::cur); // Skip filesize ...
			file.read( (char*)id, 4);  
			if( memcmp( id, "WAVE", 4))
				throw std::runtime_error( "wavfile_t::wavfile_t(): no 'WAVE' id in file '" + filename + "'.");

			unsigned int i, r, a;
			signed short c, s, f;

			// Find where fmt chunk is
			while( !file.rdstate()){
				file.read( (char*)id, 4);
				if( !memcmp( id, "fmt ", 4))
					goto fmt_code;  // Ouch a goto!
				file.read( (char*)&i, sizeof( unsigned int));
				file.seekg( i, std::ios::cur);
			}
			throw std::runtime_error( "wavfile_t::wavfile_t(): no 'fmt ' chunk in file '" + filename + "'.");

		fmt_code:
			// Get important parameters
			file.read( (char*)&i, sizeof( unsigned int));
			file.read( (char*)&f, sizeof( signed short));
			file.read( (char*)&c, sizeof( signed short));
			file.read( (char*)&r, sizeof( unsigned int));
			file.read( (char*)&a, sizeof( unsigned int));
			file.seekg( 2, std::ios::cur); // Skip junk parameter
			file.read( (char*)&s, sizeof( signed short));
			channels = c;
			samplerate = r;
			bits = s;

			// Get sample format
			if( (f > 1) && (s != 16 && s != 8 && s != 32))
				throw std::runtime_error( "wavfile_t::wavfile_t(): unexpected format (not PCM float, 16-bit or 8-bit) in file '" + filename + "'.");

			// Move to first sample
			file.seekg( i-16, std::ios::cur);
			while( !file.rdstate()){
				file.read( (char*)id, 4);
				if( !memcmp( id, "data", 4))
					goto data_code;  // Ouch a goto!
				file.read( (char*)&i, sizeof( unsigned int));
				file.seekg( i, std::ios::cur);
			}
			throw std::runtime_error( "wavfile_t::wavfile_t(): no 'data' chunk in file '" + filename + "'.");

		data_code:
			// Read number of samples
			file.read( (char*)&frames, sizeof( unsigned int));
			frames /= (bits/8) / channels;

			file_type = 1;

		// Open an output file
		}else if( mode == std::ios::out){
			file.open( filename.c_str(), std::ios::out | std::ios::binary);
			if (!file)
				throw std::runtime_error( "wavfile_t::wavfile_t(): problem writing to file '" + filename + "'.");
			unsigned int i, s;

			// Write header
			file.write( "RIFF", 4);
			i = 0; file.write( (char*)&i, 4); // file size, write in the end
			file.write( "WAVE", 4);
			file.write( "fmt ", 4);
			i = 16; file.write( (char*)&i, 4); // chunk length
			s = 1; file.write( (char*)&s, 2); // PCM 
			file.write( (char*)&channels, 2); // channels
			i = (unsigned int)samplerate; file.write( (char*)&i, 4); // sample rate
			i = (unsigned int)samplerate*2*channels; file.write( (char*)&i, 4); // bytes per second
			s = 2*channels; file.write( (char*)&s, 2); // bytes per frame
			s = 16; file.write( (char*)&s, 2); // bits per sample

			// Start with the data chunk
			file.write( "data", 4);
			i = 0; file.write( (char*)&i, 4); // number of sample bytes, write in the end

			file_type = 2;
		}
	}

	~wavfile_t()
	{
		// Write sizes
		if( file_type == 2){
			file.seekp( 0, std::ios::end);
			unsigned int i = file.tellp();
			i -= 4;
			file.seekp( 4, std::ios::beg);
			file.write( (char*)&i, 4);
			file.seekp( 40, std::ios::beg);
			i -= 8+24+8;
			file.write( (char*)&i, 4);
		}
	}

	// Read a mono buffer
	template <class T>
	void read_mono( array<T> &p)
	{
		if( file_type == 2)
			throw std::runtime_error( "wavfile_t::read_mono(): file opened write-only");

		// Read only the first channel
		if( bits == 16){
			short s[16];
			for( size_t i = 0 ; (i < p.size()) && !file.eof() ; ++i){
				file.read( (char*)s, channels*sizeof( short));
				p(i) = double( s[0])/32768.0;
			}
		}else if( bits == 8){
			signed char s[16];
			for( size_t i = 0 ; (i < p.size()) && !file.eof() ; ++i){
				file.read( (char*)s, channels*sizeof( char));
				p(i) = double( s[0])/128.0;
			}
		}else if( bits == 32){
			float s[16];
			for( size_t i = 0 ; (i < p.size()) && !file.eof() ; ++i){
				file.read( (char*)s, channels*sizeof( float));
				p(i) = double( s[0]);
			}
		}
	}

	// Write a (mono?) buffer
	template <class T>
	void write( const array<T> &p)
	{
		for( size_t i = 0 ; i < p.size() ; ++i){
			const short s = 32767*p(i);
			file.write( (char*)&s, sizeof( short));
		}
	}
};
