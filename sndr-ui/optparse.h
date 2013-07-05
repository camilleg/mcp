// optparse.h -- Functions to parse option string
// language: C++
// author  : Paris Smaragdis
// copyright: (c) 2009 Adobe Systems Incorporated, All Rights Reserved

#ifndef __OPTPARSE_H
#define __OPTPARSE_H

#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>  // for strncmp() and strcmp()
#include <cstdlib>
#include <vector>
#include <typeinfo> // for typeid() with gcc 4.4+

//
// Option parsing routine
//

template <class T> T optassign( const char *o) { std::cout << typeid( T).name() << std::endl; throw std::runtime_error( "Unsupported data type in optassign"); }
template <> int optassign<int>( const char *o) { return atoi( o); }
template <> double optassign<double>( const char *o) { return atof( o); }
template <> std::string optassign<std::string>( const char *o) { return std::string( o); }

template <class T>
T getoption( std::string opt, int argc, const char **argv, T iv = T(), const char *descr = NULL)
{
	T v = iv;
	for( int i = 0 ; i < argc ; i++){
		// Find option flag
		if( !strncmp( argv[i], opt.c_str(), opt.size())){

			// Assign accordingly
			if( strcmp( argv[i], opt.c_str()) == 0)
				v = optassign<T>( argv[i+1]);
			else
				v = optassign<T>( argv[i]+opt.size());

			if( descr)
				std::cout << descr << " is " << v << std::endl;
			else
				std::cout << "Flag " << opt << " is " << v << std::endl;
		}
	}
	return v;
}

template <class T>
array<T> mgetoption( std::string opt, int argc, const char **argv, const char *descr = NULL)
{
	array<T> v;
	for( int i = 0 ; i < argc ; i++){
		// Find option flag
		if( !strncmp( argv[i], opt.c_str(), opt.size())){
			
			// Assign accordingly
			if( strcmp( argv[i], opt.c_str()) != 0)
				v.push_back( optassign<T>( argv[i]+opt.size()));
			for( int j = i+1 ; j < argc & argv[j][0] != '-' ; j++)
				v.push_back( optassign<T>( argv[j]));

			if( descr)
				std::cout << descr << " is ";
			else
				std::cout << "Flag " << opt << " is ";
			for( int j = 0 ; j < v.size() ; j++)
				std::cout << v(j) << ", ";
			std::cout << std::endl;
		}
	}

	// By default leave empty
	return v;
}

#endif
