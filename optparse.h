// optparse.h -- Functions to parse option string
// language: C++
// author  : Paris Smaragdis

#ifndef __OPTPARSE_H
#define __OPTPARSE_H

#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>  // for strncmp() and strcmp()
#include <cstdlib>
#include <vector>
#include <typeinfo> // for typeid() with gcc 4.4+

// Option parser

template <class T> T    optassign             ( const char * )
{
	std::cout << typeid( T).name() << std::endl;
	throw std::runtime_error( "optassign: unsupported data type.");
}
template <> int         optassign<int        >( const char *o) { return atoi( o);        }
template <> double      optassign<double     >( const char *o) { return atof( o);        }
template <> std::string optassign<std::string>( const char *o) { return std::string( o); }

// Get argv after "opt".
// If multiple opt's, use the last one.
// Likely bug if the last one isn't followed by a value ("a.out -aFlag a -bFlag b -cFlag").
template <class T>
T getoption( const std::string& opt, int argc, const char **argv, T iv = T(), const char *descr = NULL)
{
	T v = iv;
	for( int i = 0 ; i < argc ; ++i){
		// Find option flag
		if( !strncmp( argv[i], opt.c_str(), opt.size())){
			// Assign accordingly
			v = optassign<T>( strcmp( argv[i], opt.c_str()) == 0 ? argv[++i] : argv[i]+opt.size());

			if( descr)
				std::cout << descr << " is " << v << std::endl;
			else
				std::cout << "Flag " << opt << " is " << v << std::endl;
		}
	}
	return v;
}

// Accumulate argvs after "opt" until the next argv that begins with a '-'.
// Unlike valarray<T> or vector<T>, array<T> uses contiguous memory slots
// (even though that doesn't apply when T is std::string).
template <class T>
array<T> mgetoption( const std::string& opt, int argc, const char **argv, const char *descr = NULL)
{
	array<T> v;
	for( int i = 0 ; i < argc ; ++i){
		// Find option flag
		if( !strncmp( argv[i], opt.c_str(), opt.size())){
			// Assign accordingly
			if( strcmp( argv[i], opt.c_str()) != 0)
				v.push_back( optassign<T>( argv[i]+opt.size()));
			for( int j = ++i ; j < argc && argv[j][0] != '-' ; ++j) {
				v.push_back( optassign<T>( argv[j]));
				++i;
			}

			if( descr)
				std::cout << descr << " is ";
			else
				std::cout << "Flag " << opt << " is ";
			for( size_t j = 0 ; j < v.size() ; ++j)
				std::cout << v(j) << ((j < v.size()-1) ? ", " : ".\n");
		}
	}

	// By default leave empty
	return v;
}

#endif
