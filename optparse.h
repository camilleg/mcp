// optparse.h -- Functions to parse option string
// language: C++
// author  : Paris Smaragdis

#ifndef __OPTPARSE_H
#define __OPTPARSE_H

#include <iostream>
#include <stdexcept>
#include <cstring>  // for strncmp() and strcmp()
#include <cstdlib>
#include <typeinfo> // for typeid() with gcc 4.4+

template <class T> T    optassign             ( const char * )
{
	std::cout << typeid( T).name() << std::endl;
	throw std::runtime_error( "optassign: unsupported data type.");
}
template <> int         optassign<int        >( const char *o) { return atoi( o);        }
template <> double      optassign<double     >( const char *o) { return atof( o);        }
template <> std::string optassign<std::string>( const char *o) { return std::string( o); }

// Get argv after "opt", either -x42 or -x 42.
// If multiple opt's, use the last one.
// Likely bug if the last one isn't followed by a value ("a.out -aFlag a -bFlag b -cFlag").
template <class T>
T getoption( const std::string& opt, int argc, const char **argv, T iv = T(), const char *descr = NULL)
{
	T v = iv;
	for( int i = 1 ; i < argc ; ++i){
		// Find option flag
		if( !strncmp( argv[i], opt.c_str(), opt.size())){
			// Assign accordingly
			v = optassign<T>( strcmp( argv[i], opt.c_str()) == 0 ? argv[++i] : argv[i]+opt.size());

			if( descr)
				std::cout << argv[0] << ": " << descr << " is " << v << "." << std::endl;
			else
				std::cout << argv[0] << ": flag " << opt << " is " << v << "." << std::endl;
			// Don't break.  Keep looking for later flags.
		}
	}
	return v;
}

// Accumulate argvs after "opt" until the next argv that begins with a '-'.
// (This *could* return std::valarray, which also guarantees contiguous memory (C++03 26.3.2.3/3),
// or even std::vector, also contiguous (C++03 23.2.4.1).  But that would propagate a mishmash
// of array<T> and std::valarray<T>'s into the rest of the code.)
template <class T>
array<T> mgetoption( const std::string& opt, int argc, const char **argv, const char *descr = NULL)
{
	array<T> v; // empty by default
	for( int i = 1 ; i < argc ; ++i){
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
				std::cout << argv[0] << ": " << descr << " is ";
			else
				std::cout << argv[0] << ": flag " << opt << " is ";
			for( size_t j = 0 ; j < v.size() ; ++j)
				std::cout << v(j) << ((j < v.size()-1) ? ", " : ".\n");
		}
	}
	return v;
}

#endif
